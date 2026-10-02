"""Independent synthetic-data arithmetic oracle. Does NOT execute or validate R."""
import json
import re


def cigar_intervals(pos1, cigar):
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    assert "".join(n + op for n, op in tokens) == cigar
    p = pos1 - 1
    blocks, introns = [], []
    for length, op in tokens:
        n = int(length)
        if op in "M=X":
            blocks.append((p, p + n))
        if op == "N":
            introns.append((p, p + n))
        if op in "M=XDN":
            p += n
    return blocks, introns


records = [
    ("exact01", 0, 181, 60, "20M100N30M"),
    ("exact02", 0, 171, 60, "30M100N20M"),
    ("near01", 0, 183, 60, "20M100N30M"),
    ("plain01", 0, 351, 60, "50M"),
    ("lowmapq01", 0, 181, 5, "20M100N30M"),
    ("secondary01", 256, 181, 60, "20M100N30M"),
    ("supplementary01", 2048, 181, 60, "20M100N30M"),
    ("span_only01", 0, 1, 60, "50M500N50M"),
]
roi0 = (100, 450)
accepted, events = set(), []
depth = {p: 0 for p in range(*roi0)}
for name, flag, pos1, mq, cigar in records:
    if flag & 2820 or mq < 20:
        continue
    blocks, introns = cigar_intervals(pos1, cigar)
    if not any(a < roi0[1] and b > roi0[0] for a, b in blocks):
        continue
    accepted.add(name)
    for a, b in blocks:
        for p in range(max(a, roi0[0]), min(b, roi0[1])):
            depth[p] += 1
    for a, b in introns:
        events.append((name, a, b))
exact = {name for name, a, b in events if (a, b) == (200, 300)}
near = {name for name, a, b in events
        if abs(a-200) <= 5 and abs(b-300) <= 5 and (a, b) != (200, 300)}
answer = dict(denominator_reads=len(accepted), junction_reads=len({x[0] for x in events}),
              junction_exact_count=len(exact), junction_near_count=len(near),
              junction_per_100_reads=100*len(exact)/len(accepted),
              mean_read_depth=sum(depth.values())/len(depth))
assert len(accepted) == 4 and len(exact) == 2 and len(near) == 1
assert sum(depth.values()) == 200
assert not cigar_intervals(181, "20M100D30M")[1]
assert len(cigar_intervals(181, "20M100N30M100N20M")[1]) == 2
assert cigar_intervals(181, "5S10=2X8M3I100N20M2D5M")[1] == [(200, 300)]
report = {"synthetic_expected_values": answer,
          "independent_python_oracle": "PASS",
          "scope": "Synthetic arithmetic only; R, native tools and browser checks are separate."}
print(json.dumps(report, indent=2, ensure_ascii=False))
