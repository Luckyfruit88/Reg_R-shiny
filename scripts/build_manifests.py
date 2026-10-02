#!/usr/bin/env python3
"""Build Reg_Shiny viewer manifests on SCC; individual rows never go to stdout."""
import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args):
    return subprocess.check_output(args, universal_newlines=True)


def read_csv(path, required):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError("Missing columns in " + str(path))
        return list(reader)


def write_tsv(path, columns, rows):
    temp = path.with_name(path.name + ".partial")
    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(str(temp), str(path))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rna-map", required=True)
    ap.add_argument("--bam-map", required=True)
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--build", default="GRCh38", choices=["GRCh38"])
    args = ap.parse_args()
    # Private until the deploying process deliberately assigns group read access.
    os.umask(0o077)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    targets = [out / n for n in ("sample_manifest.tsv", "vcf_manifest.tsv", "manifest_build_audit.json")]
    if any(p.exists() for p in targets):
        raise ValueError("Output already exists; choose a fresh output directory")
    rna = read_csv(args.rna_map, ["NWGC_ID", "framid", "Batch_ID"])
    bams = read_csv(args.bam_map, ["NWGC_ID", "directory"])
    if len({r["NWGC_ID"] for r in rna}) != len(rna):
        raise ValueError("Duplicate RNA sample IDs require source-level resolution")
    bam_map = collections.defaultdict(set)
    for row in bams:
        bam_map[row["NWGC_ID"]].add(row["directory"])
    # Validate each chromosome header; no full VCF scans or hashes of multi-GB inputs.
    vcfs = []
    header_audit = []
    expected_ids = None
    expected_lengths = {"chr1": 248956422, "chr17": 83257441,
                        "chr21": 46709983, "chrX": 156040895}
    for chrom in ["chr" + str(n) for n in range(1, 23)] + ["chrX"]:
        path = Path(args.vcf_dir) / (chrom + ".vcf.gz")
        index = next((p for p in [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
                      if p.is_file() and os.access(str(p), os.R_OK)), None)
        if not path.is_file() or not os.access(str(path), os.R_OK) or index is None:
            raise ValueError("VCF or readable index missing: " + str(path))
        header = run(["bcftools", "view", "--no-version", "-h", str(path)])
        sample_line = next((l for l in header.splitlines() if l.startswith("#CHROM\t")), None)
        if sample_line is None:
            raise ValueError("Missing VCF sample header")
        ids = sample_line.split("\t")[9:]
        if len(set(ids)) != len(ids):
            raise ValueError("Duplicate VCF sample IDs")
        if expected_ids is None:
            expected_ids = ids
        elif expected_ids != ids:
            raise ValueError("VCF sample identities/order differ at " + chrom)
        match = re.search(r"^##contig=<ID=" + re.escape(chrom) + r",length=(\d+)(?:,|>)", header, re.M)
        if not match:
            raise ValueError("Missing exact chromosome length in VCF header")
        length = int(match.group(1))
        if chrom in expected_lengths and length != expected_lengths[chrom]:
            raise ValueError("GRCh38 sentinel contig length mismatch: " + chrom)
        vcfs.append({"chrom": chrom, "vcf": str(path.resolve()), "build": args.build,
                     "index": str(index.resolve()), "length_bp": length})
        stat = path.stat()
        header_audit.append({"chrom": chrom, "vcf": str(path), "size_bytes": stat.st_size,
                             "mtime_ns": stat.st_mtime_ns, "index_sha256": sha(index),
                             "header_sha256": hashlib.sha256(header.encode()).hexdigest(),
                             "sample_count": len(ids), "length_bp": length})
    valid_ids = set(expected_ids)
    samples = []
    exclusions = collections.Counter()
    excluded_rows = []
    for row in rna:
        sid, fid = row["NWGC_ID"], row["framid"]
        reason = None
        if not fid or fid not in valid_ids:
            reason = "NO_VCF_MATCH"
        else:
            candidates = [Path(d) / (sid + ".accepted_hits.merged.markeddups.recal.bam")
                          for d in bam_map[sid]]
            candidates = list({str(p.resolve()): p for p in candidates if p.is_file()}.values())
            if len(candidates) != 1:
                reason = "NO_BAM" if not candidates else "AMBIGUOUS_BAM"
            else:
                bam = candidates[0]
                indexes = [Path(str(bam) + ".bai"), bam.with_suffix(".bai"),
                           Path(str(bam) + ".csi"), bam.with_suffix(".csi")]
                if not os.access(str(bam), os.R_OK):
                    reason = "BAM_NOT_READABLE"
                elif not any(p.is_file() and os.access(str(p), os.R_OK) for p in indexes):
                    reason = "NO_READABLE_BAM_INDEX"
        if reason:
            exclusions[reason] += 1
            excluded_rows.append({"sample_id": sid, "vcf_sample": fid, "reason": reason})
        else:
            samples.append({"vcf_sample": fid, "bam": str(bam.resolve()), "build": args.build,
                            "sample_id": sid, "batch": row["Batch_ID"]})
    if not samples or len({r["vcf_sample"] for r in samples}) != len(samples):
        raise ValueError("No samples or repeated DNA participants: resolve before viewer deployment")
    if len({r["bam"] for r in samples}) != len(samples):
        raise ValueError("Repeated BAM identity requires source-level resolution")
    # Header-only BAM build evidence, one representative per source directory.
    representative = {}
    for row in samples:
        representative.setdefault(str(Path(row["bam"]).parent), row["bam"])
    bam_build_audit = []
    for directory, bam in sorted(representative.items()):
        header = run(["samtools", "view", "-H", bam])
        lengths = {}
        for line in header.splitlines():
            if line.startswith("@SQ\t"):
                fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
                lengths[fields["SN"]] = int(fields["LN"])
        if any(lengths.get(chrom) != length for chrom, length in expected_lengths.items()):
            raise ValueError("Representative BAM GRCh38 sentinel contig mismatch: " + directory)
        bam_build_audit.append({"bam_directory": directory, "sentinel_contig_lengths": expected_lengths,
                                "representatives_checked": 1})
    write_tsv(targets[0], ["vcf_sample", "bam", "build", "sample_id", "batch"], samples)
    write_tsv(targets[1], ["chrom", "vcf", "build", "index", "length_bp"], vcfs)
    write_tsv(out / "sample_manifest_exclusions.tsv", ["sample_id", "vcf_sample", "reason"], excluded_rows)
    audit = {"status": "PASS", "build": args.build, "rna_rows": len(rna), "bam_map_rows": len(bams),
             "vcf_samples": len(expected_ids), "matched_samples": len(samples),
             "chromosomes": len(vcfs), "excluded_counts": dict(exclusions),
             "batch_counts": dict(collections.Counter(r["batch"] for r in samples)),
             "sources": [{"path": p, "sha256": sha(p)} for p in (args.rna_map, args.bam_map)],
             "vcf_headers": header_audit, "bam_build": bam_build_audit,
             "manifest_sha256": {p.name: sha(p) for p in targets[:2]},
             "scope": "All exact RNA/DNA matches with unique indexed BAM; no expression/covariate/8973 filter",
             "build_evidence_limit": "Sentinel contig lengths; BAM headers checked per source root. Runtime verifies every chosen BAM.",
             "index_validation": "FHS preparation checks index existence/readability. Runtime regional analysis opens every selected BAM index; this setup is not a full index-consistency audit.",
             "data_access": "Every user requires existing sequencing and mtdna-alcohol group authorization"}
    temp = targets[2].with_name(targets[2].name + ".partial")
    temp.write_text(json.dumps(audit, indent=2) + "\n")
    os.replace(str(temp), str(targets[2]))
    print(json.dumps({k: audit[k] for k in ("status", "build", "rna_rows", "vcf_samples", "matched_samples", "chromosomes", "excluded_counts", "batch_counts")}))


if __name__ == "__main__":
    main()
