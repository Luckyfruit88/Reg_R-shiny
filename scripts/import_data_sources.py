#!/usr/bin/env python3
"""Validate explicit custom sources and write private Reg_Shiny manifests.

Python >=3.6. No source data/index is created, changed or inferred from filenames.
Individual sample mappings are written only under the private output directory.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


def run(tool, args):
    return subprocess.check_output([tool] + list(args), universal_newlines=True, timeout=120)


def sha(path):
    value = hashlib.sha256()
    with open(str(path), "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def stat(path):
    path = Path(path).resolve()
    value = path.stat()
    return {"path": str(path), "size_bytes": value.st_size, "mtime_ns": value.st_mtime_ns}


def readable(path):
    path = Path(path).expanduser()
    if not path.is_absolute():
        raise ValueError("Source paths must be absolute.")
    if not path.is_file() or not os.access(str(path), os.R_OK):
        raise ValueError("A configured source file is unavailable or unreadable: " + str(path))
    return path.resolve()


def read_tsv(path, columns):
    with open(str(path), newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or not set(columns).issubset(reader.fieldnames):
            raise ValueError("Mapping/registry has missing or duplicate required columns.")
        rows = list(reader)
    if not rows or any(any(row.get(name) is None or not row[name].strip() for name in columns) for row in rows):
        raise ValueError("Mapping/registry must contain nonempty required values.")
    return rows


def write_tsv(path, columns, rows):
    temp = str(path) + ".partial"
    with open(temp, "w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.chmod(temp, 0o600)
    os.replace(temp, str(path))


def bam_index(bam):
    possibilities = [Path(str(bam) + ".bai"), bam.with_suffix(".bai"),
                     Path(str(bam) + ".csi"), bam.with_suffix(".csi")]
    valid = [p.resolve() for p in possibilities if p.is_file() and os.access(str(p), os.R_OK)]
    if not valid:
        raise ValueError("A selected BAM has no readable existing BAI/CSI index: " + str(bam))
    return valid[0]


def bam_header(bam, samtools):
    header = run(samtools, ["view", "-H", str(bam)])
    contigs, groups = {}, []
    for line in header.splitlines():
        if line.startswith("@SQ\t") or line.startswith("@RG\t"):
            fields = line.split("\t")
            attrs = {}
            for item in fields[1:]:
                if ":" not in item:
                    continue
                key, value = item.split(":", 1)
                if key in attrs:
                    raise ValueError("Duplicate BAM header tag.")
                attrs[key] = value
            if fields[0] == "@SQ":
                name = attrs.get("SN", "")
                if not name or name in contigs or int(attrs.get("LN", "0")) < 1:
                    raise ValueError("Invalid BAM reference dictionary.")
                contigs[name] = int(attrs["LN"])
            else:
                groups.append(attrs)
    if not contigs:
        raise ValueError("BAM has no reference dictionary.")
    return contigs, groups, hashlib.sha256(header.encode()).hexdigest()


def reject_duplicate_vcf_contigs(path):
    # bcftools may canonicalize repeated VCF header declarations. Check the
    # original text header as well, stopping at #CHROM before any variant rows.
    # BCF uses its native parsed reference dictionary instead.
    if str(path).lower().endswith(".bcf"):
        return
    with open(str(path), "rb") as source:
        compressed = source.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    names, header_bytes = set(), 0
    with opener(str(path), "rt", encoding="utf-8") as source:
        for line in source:
            header_bytes += len(line.encode("utf-8"))
            if header_bytes > 64 * 1024 * 1024:
                raise ValueError("VCF header exceeds 64 MiB.")
            if line.startswith("#CHROM"):
                return
            if line.startswith("##contig=<"):
                match = re.search(r"(?:<|,)ID=([^,>]+)", line)
                if match:
                    name = match.group(1)
                    if name in names:
                        raise ValueError("VCF has duplicate exact contig names in its original header.")
                    names.add(name)
            elif not line.startswith("#"):
                raise ValueError("VCF sample header is missing before data rows.")
    raise ValueError("VCF sample header is missing.")


def vcf_header(path, bcftools):
    if not re.search(r"\.(vcf\.(gz|bgz)|bcf)$", str(path), re.I):
        raise ValueError("Use an indexed bgzip VCF or BCF file.")
    indexes = [Path(str(path) + suffix) for suffix in (".csi", ".tbi")]
    index = next((p.resolve() for p in indexes if p.is_file() and os.access(str(p), os.R_OK)), None)
    if index is None:
        raise ValueError("A configured VCF/BCF has no readable existing index.")
    before = [stat(path), stat(index)]
    reject_duplicate_vcf_contigs(path)
    header = run(bcftools, ["view", "--no-version", "-h", str(path)])
    if before != [stat(path), stat(index)]:
        raise ValueError("VCF/index changed while its header was being inspected.")
    sample_line = [line for line in header.splitlines() if line.startswith("#CHROM\t")]
    if len(sample_line) != 1:
        raise ValueError("VCF header is missing its sample columns.")
    samples = sample_line[0].split("\t")[9:]
    if not samples or len(set(samples)) != len(samples) or any(not name or re.search(r"\s", name) for name in samples):
        raise ValueError("VCF must contain unique, nonempty exact sample identifiers.")
    contigs = {}
    for line in header.splitlines():
        if not line.startswith("##contig=<"):
            continue
        match = re.search(r"(?:<|,)ID=([^,>]+)", line)
        length = re.search(r"(?:<|,)length=(\d+)", line)
        if match:
            name = match.group(1)
            if name in contigs:
                raise ValueError("VCF has duplicate exact contig names.")
            contigs[name] = int(length.group(1)) if length else None
    if not contigs:
        raise ValueError("VCF header must declare its chromosome names.")
    return {"samples": samples, "contigs": contigs, "index": index,
            "source_states": before, "header_sha256": hashlib.sha256(header.encode()).hexdigest()}


def indexed_vcf_contigs(path, info, bcftools):
    """Discover data-bearing contigs using only existing index statistics.

    Whole-genome headers often retain unused HLA/alternate-reference names.
    Those declarations do not mean the individual chromosome VCF has records
    there. Unknown record counts are not replaced by a full-file scan.
    """
    text = run(bcftools, ["index", "-s", str(path)])
    records, seen = {}, set()
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError("Unexpected VCF index statistics; provide an explicit VCF registry.")
        chrom, length_text, count_text = fields
        if chrom in seen:
            raise ValueError("VCF index statistics contain duplicate exact contig names.")
        seen.add(chrom)
        if not re.match(r"^[0-9]+$", count_text):
            raise ValueError("VCF index record counts are unavailable; provide an explicit VCF registry.")
        count = int(count_text)
        if count == 0:
            continue
        if chrom not in info["contigs"]:
            raise ValueError("A data-bearing indexed contig is absent from the exact VCF header.")
        if not re.match(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$", chrom):
            raise ValueError("A data-bearing indexed contig has an unsupported query name; provide an explicit registry selecting supported contigs. Names are never changed.")
        header_length = info["contigs"][chrom]
        index_length = int(length_text) if re.match(r"^[0-9]+$", length_text) else None
        if header_length is not None and index_length is not None and header_length != index_length:
            raise ValueError("VCF header and index contig lengths differ.")
        records[chrom] = count
    if not records:
        raise ValueError("The VCF index has no contigs with recorded variants; provide an explicit registry if an empty declared region is intentional.")
    if info["source_states"] != [stat(path), stat(info["index"])]:
        raise ValueError("VCF/index changed during indexed contig discovery.")
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", default="")
    parser.add_argument("--vcf-registry", default="")
    parser.add_argument("--mapping", default="")
    parser.add_argument("--bam", default="")
    parser.add_argument("--bam-dir", default="")
    parser.add_argument("--mapping-mode", choices=("explicit", "read_group"), default="explicit")
    parser.add_argument("--build", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    os.umask(0o077)
    if not re.match(r"^[A-Za-z0-9_.-]+$", args.build) or args.build.lower() in ("unknown", "na", "none"):
        raise ValueError("Declare the reference assembly explicitly.")
    if args.vcf and args.vcf_registry:
        raise ValueError("Select one VCF/BCF file or a VCF registry, not both.")
    if args.bam and args.bam_dir:
        raise ValueError("Select an exact BAM or a BAM directory, not both.")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Import output directory must be empty.")
    chosen = []
    if args.bam:
        chosen = [readable(args.bam)]
    elif args.bam_dir:
        directory = Path(args.bam_dir).expanduser()
        if not directory.is_absolute() or not directory.is_dir():
            raise ValueError("BAM directory must be an existing absolute directory.")
        chosen = sorted(readable(str(p)) for p in directory.iterdir() if p.is_file() and p.suffix.lower() == ".bam")
        if not chosen:
            raise ValueError("The chosen directory contains no direct BAM files.")
    registry, cache, vcf_states = [], {}, []
    if args.vcf:
        path = readable(args.vcf)
        cache[str(path)] = vcf_header(path, args.bcftools)
        indexed = indexed_vcf_contigs(path, cache[str(path)], args.bcftools)
        cache[str(path)]["indexed_record_counts"] = indexed
        registry = [{"chrom": chrom, "vcf": str(path), "build": args.build}
                    for chrom in sorted(indexed)]
    elif args.vcf_registry:
        registry = read_tsv(readable(args.vcf_registry), ("chrom", "vcf", "build"))
        if len(set(row["chrom"] for row in registry)) != len(registry):
            raise ValueError("VCF registry must have one row per exact chromosome.")
    for row in registry:
        if row["build"] != args.build:
            raise ValueError("Registry and requested reference builds differ.")
        path = readable(row["vcf"]); row["vcf"] = str(path)
        if str(path) not in cache:
            cache[str(path)] = vcf_header(path, args.bcftools)
        info = cache[str(path)]
        if row["chrom"] not in info["contigs"]:
            raise ValueError("Registry chromosome is absent from its exact VCF header.")
        if not re.match(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$", row["chrom"]):
            raise ValueError("A selected registry contig has an unsupported query name; chromosome names are never changed.")
        run(args.bcftools, ["query", "--regions-overlap", "0", "-r", row["chrom"] + ":1-1",
                           "-f", "%CHROM\\t%POS\\n", str(path)])
        if info["source_states"] != [stat(path), stat(info["index"])]:
            raise ValueError("VCF/index changed during the selected-contig index-open probe.")
        row["length_bp"] = info["contigs"][row["chrom"]]
        row["index"] = str(info["index"])
    for path, info in sorted(cache.items()):
        vcf_states.extend(info["source_states"])
    vcf_samples = set(sample for info in cache.values() for sample in info["samples"])
    mappings, exclusions = [], []
    if registry:
        if args.mapping_mode == "explicit":
            if not args.mapping:
                raise ValueError("DNA/RNA comparison requires an explicit mapping TSV or validated read-group matching.")
            rows = read_tsv(readable(args.mapping), ("vcf_sample", "bam"))
            for row in rows:
                if row.get("build", args.build) != args.build:
                    raise ValueError("Mapping and requested reference builds differ.")
                if row["vcf_sample"] not in vcf_samples:
                    raise ValueError("An explicit mapping sample does not exactly match any VCF sample.")
                bam = readable(row["bam"])
                if chosen and bam not in chosen:
                    exclusions.append({"vcf_sample": row["vcf_sample"], "bam": str(bam), "reason": "OUTSIDE_SELECTED_BAM_SET"})
                    continue
                mappings.append({"vcf_sample": row["vcf_sample"], "bam": str(bam), "build": args.build})
            if not chosen:
                chosen = [Path(row["bam"]) for row in mappings]
        else:
            if not chosen:
                raise ValueError("Read-group matching requires an explicit BAM file or directory.")
    elif args.mapping:
        raise ValueError("A mapping TSV requires a configured DNA VCF/BCF.")
    if not chosen:
        raise ValueError("Select at least one readable indexed BAM.")
    if len(set(str(path) for path in chosen)) != len(chosen):
        raise ValueError("The selected BAM paths resolve to duplicate physical paths.")
    # Detect hard-linked aliases as well as symlinks.
    identities = [(path.stat().st_dev, path.stat().st_ino) for path in chosen]
    if len(set(identities)) != len(identities):
        raise ValueError("The selected BAM files contain duplicate hard-link identities.")
    bam_audit, bam_states = [], []
    declared_lengths = {row["chrom"]: row["length_bp"] for row in registry}
    for bam in chosen:
        if bam.suffix.lower() != ".bam":
            raise ValueError("RNA inputs must be BAM files.")
        index = bam_index(bam)
        before = [stat(bam), stat(index)]
        contigs, groups, header_hash = bam_header(bam, args.samtools)
        shared = set(contigs).intersection(declared_lengths)
        if registry and not shared:
            raise ValueError("BAM and VCF have no exact chromosome names in common; aliases are not inferred.")
        if any(declared_lengths[chrom] is not None and declared_lengths[chrom] != contigs[chrom] for chrom in shared):
            raise ValueError("BAM/VCF chromosome lengths differ; the reference assembly is incompatible.")
        # Prefer the selected DNA/RNA intersection, rather than an unrelated
        # unused header declaration. BAM-only input still requires an exact name
        # that the existing Reg_Shiny regional-query contract supports.
        probe_contigs = sorted(shared) if registry else sorted(
            chrom for chrom in contigs if re.match(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$", chrom))
        if not probe_contigs:
            raise ValueError("The BAM has no chromosome name supported by regional queries; names are never changed.")
        # Explicit regions require the existing index and cannot fall back to a
        # whole BAM scan. This proves index-open, not full index consistency.
        run(args.samtools, ["view", "-c", str(bam), probe_contigs[0] + ":1-1"])
        if registry and args.mapping_mode == "read_group":
            if not groups or any(not row.get("ID") or not row.get("SM") for row in groups) or len({row["ID"] for row in groups}) != len(groups):
                raise ValueError("Every BAM read group needs a unique ID and nonempty SM for exact matching.")
            names = {row["SM"] for row in groups}
            if len(names) != 1 or next(iter(names)) not in vcf_samples:
                raise ValueError("BAM read-group SM is mixed or does not exactly match a VCF sample; provide an explicit mapping instead.")
            mappings.append({"vcf_sample": next(iter(names)), "bam": str(bam), "build": args.build})
        if before != [stat(bam), stat(index)]:
            raise ValueError("A BAM/index changed during source validation.")
        bam_states.extend(before)
        bam_audit.append({"bam": str(bam), "index": str(index), "header_sha256": header_hash,
                          "shared_contig_count": len(shared)})
    if registry and (not mappings or len({row["vcf_sample"] for row in mappings}) != len(mappings) or
                     len({row["bam"] for row in mappings}) != len(mappings)):
        raise ValueError("DNA/RNA mapping must be nonempty and one-to-one: one selected BAM per VCF sample.")
    if any(stat(row["path"]) != row for row in vcf_states):
        raise ValueError("A VCF/index changed during import.")
    write_tsv(out / "bam_choices.tsv", ("bam",), [{"bam": str(p)} for p in sorted(chosen)])
    if registry:
        write_tsv(out / "vcf_manifest.tsv", ("chrom", "vcf", "build", "index", "length_bp"), sorted(registry, key=lambda row: row["chrom"]))
        write_tsv(out / "sample_manifest.tsv", ("vcf_sample", "bam", "build"), sorted(mappings, key=lambda row: row["vcf_sample"]))
    write_tsv(out / "mapping_exclusions.tsv", ("vcf_sample", "bam", "reason"), exclusions)
    audit = {"schema": "regshiny-source-import-v1", "status": "PASS", "build": args.build,
        "mapping_mode": args.mapping_mode if registry else "BAM_ONLY", "bam_count": len(chosen),
        "matched_samples": len(mappings), "chromosomes": len(registry), "excluded_count": len(exclusions),
        "source_states": vcf_states + bam_states, "bam_headers": bam_audit,
        "vcf_headers": [{"vcf": path, "header_sha256": info["header_sha256"], "sample_count": len(info["samples"])}
                        for path, info in sorted(cache.items())],
        "registry_selection": "Single-file import uses data-bearing contigs reported by the existing index; explicit registries select exact user-specified header contigs.",
        "indexed_record_counts": {path: info["indexed_record_counts"] for path, info in sorted(cache.items())
                                  if "indexed_record_counts" in info},
        "input_files": [{"path": str(readable(path)), "sha256": sha(readable(path))}
                        for path in (args.mapping, args.vcf_registry) if path],
        "identity_policy": "Explicit mapping or exact BAM RG SM / VCF sample equality. No filename/participant-ID inference.",
        "index_validation": "Existing indexes opened through explicit one-base region queries; no full-scan fallback. This is not a proof of full index/data consistency."}
    audit_path = out / "source_import_audit.json"
    audit_path.write_text(json.dumps(audit, sort_keys=True, indent=2) + "\n")
    os.chmod(str(audit_path), 0o600)
    print(json.dumps({key: audit[key] for key in ("status", "build", "mapping_mode", "bam_count", "matched_samples", "chromosomes", "excluded_count")}))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("Source import failed: " + str(error), file=__import__("sys").stderr)
        __import__("sys").exit(2)
