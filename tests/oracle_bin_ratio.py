"""Independent arithmetic oracle for splice-site bin counts. Does NOT execute or validate R.

Reads tests/fixtures/bin_ratio_reads.json, applies the documented rules from
scratch and checks the hand-derived expected [intron, exon] counts per site:
a retained read overlaps the analysis interval with an aligned block; it counts
toward a bin when an M/=/X block overlaps the bin by >= 1 bp (N and D never
count); RF/FR keep reads inferred on the site strand (unknown never counts);
XS counts all strands.
"""
import json
import pathlib
import re

FIXTURE = pathlib.Path(__file__).with_name("fixtures") / "bin_ratio_reads.json"


def aligned_blocks(pos1, cigar):
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    assert "".join(n + op for n, op in tokens) == cigar
    ref, blocks = pos1, []          # 1-based reference position of the next base
    for length, op in tokens:
        n = int(length)
        if op in "M=X":
            blocks.append((ref, ref + n - 1))
        if op in "M=XDN":
            ref += n
    return blocks


def library_strand(flag, mode):
    """RegTools flag convention for paired/single reads; '?' if mates disagree."""
    reverse, mate_reverse = bool(flag & 16), bool(flag & 32)
    read1, read2 = bool(flag & 64), bool(flag & 128)
    flip = mode == "RF"
    a = (flip != read1) != reverse
    b = (flip != read2) != mate_reverse
    if a != b:
        return "?"
    return "+" if a else "-"


def bins(site, intron_width, exon_width):
    b = site["position"]
    forward = (site["type"] == "donor") == (site["strand"] == "+")
    if forward:
        return (b, b + intron_width - 1), (b - exon_width, b - 1)
    return (b - intron_width + 1, b), (b + 1, b + exon_width)


def overlaps(blocks, lo, hi):
    return any(s <= hi and e >= lo for s, e in blocks)


fx = json.loads(FIXTURE.read_text())
lo, hi = fx["interval"]
for case in fx["cases"]:
    retained = []
    for r in fx["reads"]:
        blocks = aligned_blocks(r["pos1"], r["cigar"])
        if not overlaps(blocks, lo, hi):
            continue
        nh = next((int(t[5:]) for t in r["tags"] if t.startswith("NH:i:")), None)
        if case["nh1_only"] and nh != 1:
            continue
        retained.append((blocks, library_strand(r["flag"], case["strand_mode"])))
    got = []
    for site in fx["sites"]:
        intron, exon = bins(site, case["intron_width"], case["exon_width"])
        eligible = [blk for blk, strand in retained
                    if case["strand_mode"] == "XS" or strand == site["strand"]]
        got.append([sum(overlaps(blk, *intron) for blk in eligible),
                    sum(overlaps(blk, *exon) for blk in eligible)])
    assert got == case["expected"], (case, got)
assert aligned_blocks(981, "20M500N30M") == [(981, 1000), (1501, 1530)]
assert aligned_blocks(1001, "10S5M3I5M") == [(1001, 1005), (1006, 1010)]
assert library_strand(65, "RF") == "?" and library_strand(163, "RF") == "+" and library_strand(163, "FR") == "-"
print(json.dumps({"independent_bin_ratio_oracle": "PASS", "cases": len(fx["cases"])}))
