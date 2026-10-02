#!/usr/bin/env python3
"""Build/query a GENCODE v48 index. Python >=3.6; standard library only.

GTF coordinates remain 1-based inclusive. Queries return full exon models for
every overlapping transcript. They never derive introns from clipped exons.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import re
import sqlite3
import sys
import tempfile
from pathlib import Path

SCHEMA_VERSION = "regshiny_gencode_v1"
FEATURES = {"CDS", "UTR", "start_codon", "stop_codon", "Selenocysteine"}


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source_state(path):
    stat = os.stat(str(path))
    return (stat.st_size, stat.st_mtime_ns)


def attributes(text):
    # GTF quoted values can include semicolons; do not split quoted strings.
    result = {}
    position = 0
    token = re.compile(r'\s*([^\s;]+)\s+(?:"((?:[^"\\]|\\.)*)"|([^;\s]+))\s*;?')
    while position < len(text):
        if not text[position:].strip():
            break
        match = token.match(text, position)
        if not match:
            raise ValueError("Malformed GTF attributes: " + text[:160])
        value = match.group(2) if match.group(2) is not None else match.group(3)
        value = value.replace(r'\"', '"').replace(r'\\', '\\')
        result.setdefault(match.group(1), []).append(value)
        position = match.end()
    return result


def one(attrs, key, default=""):
    values = attrs.get(key, [])
    if len(values) > 1:
        raise ValueError("Repeated identifier attribute: " + key)
    return values[0] if values else default


def create_schema(con):
    con.executescript("""
      PRAGMA journal_mode=OFF;
      PRAGMA synchronous=OFF;
      PRAGMA temp_store=MEMORY;
      CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
      CREATE TABLE genes (
        gene_id TEXT PRIMARY KEY, gene_name TEXT NOT NULL, gene_type TEXT NOT NULL,
        chrom TEXT NOT NULL, start1 INTEGER NOT NULL, end1 INTEGER NOT NULL,
        strand TEXT NOT NULL, source TEXT NOT NULL);
      CREATE TABLE transcripts (
        transcript_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, gene_name TEXT NOT NULL,
        transcript_name TEXT NOT NULL, transcript_type TEXT NOT NULL,
        transcript_support_level TEXT NOT NULL, tags TEXT NOT NULL,
        chrom TEXT NOT NULL, start1 INTEGER NOT NULL, end1 INTEGER NOT NULL,
        strand TEXT NOT NULL, source TEXT NOT NULL);
      CREATE TABLE exons (
        transcript_id TEXT NOT NULL, gene_id TEXT NOT NULL, exon_id TEXT NOT NULL,
        exon_number TEXT NOT NULL, chrom TEXT NOT NULL, start1 INTEGER NOT NULL,
        end1 INTEGER NOT NULL, strand TEXT NOT NULL,
        UNIQUE(transcript_id, start1, end1));
      CREATE TABLE features (
        transcript_id TEXT NOT NULL, gene_id TEXT NOT NULL, feature TEXT NOT NULL,
        chrom TEXT NOT NULL, start1 INTEGER NOT NULL, end1 INTEGER NOT NULL,
        strand TEXT NOT NULL, frame TEXT NOT NULL);
      CREATE TABLE introns (
        transcript_id TEXT NOT NULL, gene_id TEXT NOT NULL, chrom TEXT NOT NULL,
        start1 INTEGER NOT NULL, end1 INTEGER NOT NULL, strand TEXT NOT NULL,
        intron_number INTEGER NOT NULL);
    """)


def build(args):
    if args.release != "48" or args.build != "GRCh38":
        raise ValueError("This index contract requires GENCODE release 48 and GRCh38.")
    source = Path(args.gtf).resolve()
    output = Path(args.output).resolve()
    if output.exists() and not args.replace:
        raise ValueError("Output already exists; choose another path or explicitly pass --replace.")
    if not output.parent.is_dir():
        raise ValueError("Output parent directory must already exist.")
    before = source_state(source)
    source_hash = sha256(source)
    reference_lengths, reference_metadata = {}, None
    if args.fai:
        fai = Path(args.fai).resolve()
        fai_before = source_state(fai)
        with open(str(fai), "rt", encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 2 or not fields[0] or fields[0] in reference_lengths:
                    raise ValueError("Invalid or duplicate reference FAI contig.")
                length = int(fields[1])
                if length < 1:
                    raise ValueError("Reference FAI length must be positive.")
                reference_lengths[fields[0]] = length
        reference_metadata = {"path": str(fai), "sha256": sha256(fai),
            "description": "Administrator-supplied GRCh38 FASTA index for exact-name contig length compatibility. This does not independently authenticate the GENCODE GRCh38.p14 ALL reference; absent contigs remain unknown."}
    handle, temporary = tempfile.mkstemp(prefix=".annotation_build_", suffix=".sqlite", dir=str(output.parent))
    os.close(handle)
    con = None
    try:
        con = sqlite3.connect(temporary)
        create_schema(con)
        counts = {"genes": 0, "transcripts": 0, "exons": 0, "features": 0, "introns": 0}
        headers = []
        opener = gzip.open if str(source).endswith(".gz") else open
        with opener(str(source), "rt", encoding="utf-8") as stream:
            for number, line in enumerate(stream, 1):
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#"):
                    headers.append(line)
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError("GTF line {} does not have nine columns.".format(number))
                chrom, provider, feature, start, end, score, strand, frame, raw_attrs = fields
                if feature not in ({"gene", "transcript", "exon"} | FEATURES):
                    continue
                start, end = int(start), int(end)
                if start < 1 or end < start or strand not in ("+", "-"):
                    raise ValueError("Invalid GTF coordinates/strand at line {}.".format(number))
                if chrom in reference_lengths and end > reference_lengths[chrom]:
                    raise ValueError("GTF coordinate exceeds the supplied reference FAI contig length.")
                attrs = attributes(raw_attrs)
                gene_id = one(attrs, "gene_id")
                transcript_id = one(attrs, "transcript_id")
                if not gene_id or (feature != "gene" and not transcript_id):
                    raise ValueError("Missing GTF gene/transcript ID at line {}.".format(number))
                if feature == "gene":
                    con.execute("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?)", (
                        gene_id, one(attrs, "gene_name"), one(attrs, "gene_type"), chrom,
                        start, end, strand, provider))
                    counts["genes"] += 1
                elif feature == "transcript":
                    con.execute("INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (
                        transcript_id, gene_id, one(attrs, "gene_name"), one(attrs, "transcript_name"),
                        one(attrs, "transcript_type"), one(attrs, "transcript_support_level"),
                        ";".join(attrs.get("tag", [])), chrom, start, end, strand, provider))
                    counts["transcripts"] += 1
                elif feature == "exon":
                    con.execute("INSERT INTO exons VALUES (?,?,?,?,?,?,?,?)", (
                        transcript_id, gene_id, one(attrs, "exon_id"), one(attrs, "exon_number"),
                        chrom, start, end, strand))
                    counts["exons"] += 1
                else:
                    con.execute("INSERT INTO features VALUES (?,?,?,?,?,?,?,?)", (
                        transcript_id, gene_id, feature, chrom, start, end, strand, frame))
                    counts["features"] += 1
                if number % 100000 == 0:
                    con.commit()
        description = " ".join(headers)
        if not re.search(r"\bGRCh38\b", description) or not re.search(r"\bversion\s+48\b", description, re.I):
            raise ValueError("GTF header must explicitly identify GRCh38 and version 48.")
        if not counts["genes"] or not counts["transcripts"] or not counts["exons"]:
            raise ValueError("The complete GENCODE gene/transcript/exon model is required.")
        con.executescript("""
          CREATE INDEX genes_region ON genes(chrom,start1,end1);
          CREATE INDEX transcripts_region ON transcripts(chrom,start1,end1);
          CREATE INDEX transcripts_gene ON transcripts(gene_id);
          CREATE INDEX exons_transcript ON exons(transcript_id,start1,end1);
          CREATE INDEX features_transcript ON features(transcript_id,start1,end1);
        """)
        bad = con.execute("""SELECT t.transcript_id FROM transcripts t LEFT JOIN genes g USING(gene_id)
            WHERE g.gene_id IS NULL OR t.chrom!=g.chrom OR t.strand!=g.strand
            OR t.start1<g.start1 OR t.end1>g.end1 LIMIT 1""").fetchone()
        if bad:
            raise ValueError("Transcript lacks a consistent complete gene model: " + bad[0])
        for table in ("exons", "features"):
            bad = con.execute("""SELECT e.transcript_id FROM {} e LEFT JOIN transcripts t USING(transcript_id)
              WHERE t.transcript_id IS NULL OR e.gene_id!=t.gene_id OR e.chrom!=t.chrom OR e.strand!=t.strand
              OR e.start1<t.start1 OR e.end1>t.end1 LIMIT 1""".format(table)).fetchone()
            if bad:
                raise ValueError("Feature lacks a consistent complete transcript model: " + bad[0])
        bad = con.execute("""SELECT t.transcript_id FROM transcripts t
          WHERE NOT EXISTS (SELECT 1 FROM exons e WHERE e.transcript_id=t.transcript_id) LIMIT 1""").fetchone()
        if bad:
            raise ValueError("Transcript has no exon model: " + bad[0])
        # One transcript at a time: original full genomic exon boundaries.
        rows = con.execute("""SELECT e.transcript_id,e.gene_id,e.chrom,e.start1,e.end1,e.strand,
          e.exon_number,t.start1,t.end1 FROM exons e JOIN transcripts t USING(transcript_id)
          ORDER BY e.transcript_id,e.start1,e.end1""")
        previous_id, model = None, []

        def add_introns(exons):
            try:
                numbers = [int(exon[6]) for exon in exons]
            except ValueError:
                raise ValueError("GENCODE exon_number is missing or invalid in transcript " + exons[0][0])
            expected = list(range(1, len(exons) + 1))
            if exons[0][5] == "-":
                expected.reverse()
            if numbers != expected:
                raise ValueError("Incomplete or inconsistent exon_number sequence in transcript " + exons[0][0])
            if exons[0][3] != exons[0][7] or exons[-1][4] != exons[0][8]:
                raise ValueError("Exon model does not cover full transcript boundaries: " + exons[0][0])
            pending = []
            for index in range(1, len(exons)):
                left, right = exons[index - 1], exons[index]
                if right[3] <= left[4]:
                    raise ValueError("Overlapping exon boundaries in transcript " + left[0])
                start, end = left[4] + 1, right[3] - 1
                if end >= start:
                    rank = index if left[5] == "+" else len(exons) - index
                    pending.append((left[0], left[1], left[2], start, end, left[5], rank))
            con.executemany("INSERT INTO introns VALUES (?,?,?,?,?,?,?)", pending)
            counts["introns"] += len(pending)

        for row in rows:
            if previous_id is not None and row[0] != previous_id:
                add_introns(model)
                model = []
            model.append(row)
            previous_id = row[0]
        if model:
            add_introns(model)
        con.executescript("""
          CREATE INDEX introns_transcript ON introns(transcript_id,start1,end1);
          CREATE INDEX introns_match ON introns(chrom,start1,end1,strand);
        """)
        if source_state(source) != before:
            raise ValueError("Source GTF changed while the index was being built.")
        if args.fai and source_state(fai) != fai_before:
            raise ValueError("Reference FAI changed while the index was being built.")
        metadata = {
            "schema_version": SCHEMA_VERSION, "release": "48", "build": "GRCh38", "assembly": "GRCh38.p14",
            "coordinate_system": "1-based inclusive; source GTF boundaries unchanged",
            "source_path": str(source), "source_filename": source.name,
            "source_url": args.source_url, "source_sha256": source_hash,
            "source_size_bytes": before[0], "source_header": headers,
            "reference_contig_lengths": reference_lengths, "reference_fai": reference_metadata,
            "builder_sha256": sha256(Path(__file__).resolve()), "counts": counts,
            "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "python_version": sys.version.split()[0], "sqlite_version": sqlite3.sqlite_version,
            "intron_definition": "Consecutive full-model exons in genomic order: previous exon end+1 through next exon start-1; original strand retained."
        }
        con.executemany("INSERT INTO metadata VALUES (?,?)", [(key, json.dumps(value, sort_keys=True)) for key, value in sorted(metadata.items())])
        con.commit()
        if con.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("SQLite integrity check failed.")
        con.close()
        con = None
        # Public annotation derivative; the containing application controls its
        # access policy. No source GTF is changed by indexing.
        os.chmod(temporary, 0o640)
        os.replace(temporary, str(output))
        return {"database": str(output), "metadata": metadata, "database_sha256": sha256(output)}
    finally:
        if con is not None:
            con.close()
        if os.path.exists(temporary):
            os.unlink(temporary)


def open_readonly(path):
    path = Path(path).resolve()
    con = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only=ON")
    metadata = {row["key"]: json.loads(row["value"]) for row in con.execute("SELECT key,value FROM metadata")}
    if metadata.get("schema_version") != SCHEMA_VERSION or metadata.get("release") != "48" or metadata.get("build") != "GRCh38":
        con.close()
        raise ValueError("Annotation index schema/release/build does not match GENCODE 48 / GRCh38.")
    return con, metadata


def metadata_only(args):
    before = source_state(args.db)
    con, metadata = open_readonly(args.db)
    try:
        metadata["contigs"] = [row[0] for row in con.execute("SELECT DISTINCT chrom FROM genes ORDER BY chrom")]
    finally:
        con.close()
    if source_state(args.db) != before:
        raise ValueError("Annotation index changed while being read.")
    return metadata


def query(args):
    if not re.match(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$", args.chrom):
        raise ValueError("Invalid chromosome name.")
    if args.start < 1 or args.end < args.start or args.end - args.start + 1 > 250000:
        raise ValueError("The annotation interval must be 1-based inclusive and at most 250,000 bp.")
    if args.max_transcripts < 1 or args.max_records < 1:
        raise ValueError("Annotation record limits must be positive.")
    before = source_state(args.db)
    con, metadata = open_readonly(args.db)
    try:
        con.execute("BEGIN")
        if args.build != metadata["build"]:
            raise ValueError("Variant and annotation reference assemblies differ; no liftover is performed.")
        if con.execute("SELECT 1 FROM genes WHERE chrom=? LIMIT 1", (args.chrom,)).fetchone() is None:
            raise ValueError("Chromosome is not present in this annotation; no automatic name conversion is performed.")
        region = (args.chrom, args.end, args.start)
        count = con.execute("SELECT count(*) FROM transcripts WHERE chrom=? AND start1<=? AND end1>=?", region).fetchone()[0]
        if count > args.max_transcripts:
            raise ValueError("Annotation query exceeds the transcript limit ({} > {}); select a smaller interval. No models were truncated.".format(count, args.max_transcripts))
        tx = [dict(row) for row in con.execute("SELECT * FROM transcripts WHERE chrom=? AND start1<=? AND end1>=? ORDER BY transcript_id", region)]
        genes = [dict(row) for row in con.execute("SELECT * FROM genes WHERE chrom=? AND start1<=? AND end1>=? ORDER BY gene_id", region)]
        answer = {"genes": genes, "transcripts": tx, "exons": [], "features": [], "introns": []}
        used = len(genes) + len(tx)
        if used > args.max_records:
            raise ValueError("Annotation query exceeds the record limit; select a smaller interval. No biological tables were truncated.")
        ids = [row["transcript_id"] for row in tx]
        for table in ("exons", "features", "introns"):
            for start in range(0, len(ids), 400):
                batch = ids[start:start + 400]
                placeholders = ",".join("?" for _ in batch)
                sql = "SELECT * FROM {} WHERE transcript_id IN ({}) ORDER BY transcript_id,start1,end1".format(table, placeholders)
                for row in con.execute(sql, batch):
                    used += 1
                    if used > args.max_records:
                        raise ValueError("Annotation query exceeds the {}-record limit; select a smaller interval. No biological tables were truncated.".format(args.max_records))
                    answer[table].append(dict(row))
        for row in answer["introns"]:
            row["intron_start1"], row["intron_end1"] = row["start1"], row["end1"]
            row["donor1"] = row["start1"] if row["strand"] == "+" else row["end1"]
            row["acceptor1"] = row["end1"] if row["strand"] == "+" else row["start1"]
        for table in ("exons", "features", "introns"):
            answer[table].sort(key=lambda row: (row["transcript_id"], row["start1"], row["end1"], row.get("feature", "")))
        metadata["query"] = {"chrom": args.chrom, "start1": args.start, "end1": args.end,
                             "full_transcript_models": True, "truncated": False,
                             "max_transcripts": args.max_transcripts, "max_records": args.max_records}
        metadata["returned_counts"] = {table: len(answer[table]) for table in answer}
        answer["metadata"] = metadata
    finally:
        con.close()
    if source_state(args.db) != before:
        raise ValueError("Annotation index changed while being read.")
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command")
    b = subs.add_parser("build")
    b.add_argument("--gtf", required=True)
    b.add_argument("--output", required=True)
    b.add_argument("--release", default="48")
    b.add_argument("--build", default="GRCh38")
    b.add_argument("--source-url", default="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/")
    b.add_argument("--fai", help="Optional administrator-supplied GRCh38 FASTA index; exact-name lengths only.")
    b.add_argument("--replace", action="store_true")
    m = subs.add_parser("metadata")
    m.add_argument("--db", required=True)
    q = subs.add_parser("query")
    q.add_argument("--db", required=True)
    q.add_argument("--chrom", required=True)
    q.add_argument("--start", type=int, required=True)
    q.add_argument("--end", type=int, required=True)
    q.add_argument("--build", default="GRCh38")
    q.add_argument("--max-transcripts", type=int, default=5000)
    q.add_argument("--max-records", type=int, default=200000)
    args = parser.parse_args()
    if not args.command:
        parser.error("Choose build, metadata or query.")
    try:
        value = {"build": build, "metadata": metadata_only, "query": query}[args.command](args)
        print(json.dumps(value, sort_keys=True, separators=(",", ":")))
    except Exception as error:
        print("Annotation error: " + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
