# Metric definitions

## Coordinates and counting units

The interface uses 1-based, inclusive skipped-intron coordinates. For POS=181 and CIGAR=20M100N30M, the intron is 201-300; its internal BED interval is [200,300). RegTools BED12 outer coordinates include anchors. Intron start0 = chromStart + first blockSize; intron end0 = chromEnd - second blockSize.

Reads are retained primary alignment records. Paired ends count separately, not as fragments or UMIs. Summary metrics combine strands, including unknown strand; junction tables retain strand information.

## Metrics

| Metric | Definition |
|---|---|
| Denominator reads | Retained alignments with an M, = or X block overlapping the fixed analysis interval; N-only and D-only overlaps are excluded. |
| Junction reads | Union of read IDs supporting at least one accepted RegTools junction inside the interval. |
| Exact count | Read IDs supporting the exact two target boundaries, combined across strands. |
| Near count | Union of read IDs supporting non-exact junctions whose two boundaries are each within the selected tolerance. |
| Junctions per 100 reads | 100 × exact count / denominator reads; NA for a zero denominator. This is not PSI. |
| Mean read depth | Mean samtools depth over all positions, including introns and zero coverage. |

A multi-junction read can contribute to both exact and near counts. The exported exact_or_near_reads is their union; do not assume the two counts can be added.
The target and tolerance must leave at least one base on each side inside the loaded interval. Near matching requires both boundaries, not just a shared donor or acceptor.

## Filtering and audit

Defaults: MAPQ >= 20, base quality >= 0, anchor >= 8 bp, intron length 70-500000 bp. Always exclude UNMAP, SECONDARY, QCFAIL and SUPPLEMENTARY (mask 2820). Marked duplicates remain unless explicitly excluded. Requiring NH==1 also excludes missing NH. MAPQ=255 means unavailable mapping quality, not necessarily unique mapping.

In single-BAM analysis, base quality affects depth only. In variant comparison it also filters RNA base observations; junction counting is unchanged. M, = and X contribute to depth; N and D do not. Overlapping paired ends count twice. Flank means use up to 50 bp on either side, without GTF-based exon validation.

RegTools applies its native anchor and strandedness rules. Per-alignment CIGAR support is checked against native BED12 scores; mismatches stop the analysis. Missing XS remains unknown strand. Verify RF/FR using library-specific controls.

## Limits and outputs

Limits: 250 kb, 100000 candidate alignments and 256 MiB regional SAM per analysis. Exceeding a limit produces an error rather than subsampling. External commands have a 600-second timeout.
Target, tolerance and display-window changes reuse loaded results. Filters and analysis intervals require a new run. Display zoom does not change the denominator.

Outputs: junction_metrics.csv, junctions.csv, depth.tsv, junction_read_evidence.csv and parameters_and_audit.json. Evidence omits SEQ and QUAL but can contain QNAMEs, sample names and private paths. Keep it within the authorized research environment.

Index timestamps are warning signals only. Successful regional queries do not validate the entire index or BAM. No-junction results do not establish absence of splicing. CRAM, formal PSI, UMI counts and strand-specific denominators are not implemented. GENCODE v48 is an independent reference overlay. Genotype comparison is a descriptive extension, not a differential-splicing test.

## Genotype comparison

DNA groups come from the exact VCF CHROM/POS/REF/ALT record and the configured sample identity mapping. A different ALT or reference build is a different record. No liftover or strand/dosage flip is inferred. Phase remains in raw GT; canonical group labels sort the native allele indices. Fully or partially missing calls do not enter called-genotype RNA groups. Haploid calls are not converted into diploid genotypes.

Group depth is the arithmetic mean at each position across successfully analyzed BAMs, including successful zero-coverage samples. Junction count sums and means use those same successful samples; a successful sample without a junction contributes zero. Failed BAMs are reported as failed and excluded from these denominators, not replaced by zero. A group with no successful BAMs has unavailable evidence. Read counts are not library-size normalized, and genotype groups can differ in batch, RNA quality or expression. These plots help inspect evidence; they do not establish a causal splicing effect.

For SNVs, RNA A/C/G/T/N observations are counted from the exact retained SAM snapshot, using CIGAR M/=/X to locate the reference position. Insertions/soft clips advance the query; deletions and skipped introns provide no nucleotide observation. Reverse-strand SAM SEQ is already in reference orientation. Missing sequence/quality and low-quality bases are excluded and audited. Paired ends count separately. Group and sample audit tables retain excluded skipped/deleted bases, low base quality and missing sequence/quality. A group with no valid observations is labelled separately from positive observations. These observations are not a genotype call, RNA editing call, or proof of allele-specific splicing. For non-SNVs, allele evidence is unavailable rather than zero.

In default all-matched-BAM mode, every explicitly linked BAM with a complete DNA call is selected; no numerical sample cap applies. Optional preview selects deterministically by DNA sample identifier within each observed called genotype, capped by the displayed request. Reproducibility of a preview does not make it representative. The selected sample IDs, source metadata, parameter snapshot, native tool versions, failures and timestamps are available in the audit. Source files are read-only; runtime temporary files are private to the launching user/job.

Native semantics: [bcftools query](https://samtools.github.io/bcftools/howtos/query.html), [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), [samtools depth](https://www.htslib.org/doc/samtools-depth.html).

DNA depth or genotype-quality values are not inferred from a header declaration. The configured SCC test records contain GT only; RNA depth is computed from BAMs and must not be read as DNA DP. The viewer does not reapply DNA DP/GQ quality thresholds.

## All-matched-BAM accounting and recovery

`available_n` is the number of explicitly mapped RNA BAMs in a DNA genotype group. The mapping enforces one unique BAM per unique VCF sample ID. In the FHS preset these IDs map to unique participants; custom inputs must not be assumed to represent independent people without their own metadata. In all mode, `selected_n == available_n` for every complete called genotype; missing/partial GT groups remain selected=0. During execution, selected = completed + pending and completed = success + failure. At finalization pending must be zero.

`analyzed_n` is the number of successful BAM analyses and remains the denominator of mean depth and mean junction support. Thus a successful zero-coverage sample contributes zero; a failed BAM contributes unavailable evidence and is counted in `failed_n`. Completing every attempt is distinct from successfully measuring every sample. All mode does not turn these descriptive comparisons into association tests or correct library/batch confounding.

The per-sample checkpoint stores native metrics, audit and compact depth/junction contributions. Its checksum and frozen code/input identity are verified before reuse. A fresh aggregate is reconstructed exactly once from completed sample contributions after an interruption, preventing double counting. Complete results additionally require matching VCF/BAM/index state and terminal scheduler success. UI disconnection does not cancel the submitted scheduler job.

## GENCODE v48 overlay semantics

GENCODE features use 1-based inclusive coordinates, matching the application's displayed genomic coordinates. For adjacent exons [a,b] and [c,d] of the same complete transcript, the skipped intron is [b+1,c-1]. Negative-strand transcript direction changes donor/acceptor labelling, not the genomic bounds. Introns are derived before any display clipping. Exon numbers must be consecutive in transcript 5′→3′ order and the complete exon model must span the declared transcript; incomplete models cannot create a false annotation match.

Junction matching preserves exact chromosome, both boundaries, strand and versioned gene/transcript IDs. One observed junction can match multiple transcripts; the detailed matches retain this relationship. Unknown RNA strand is not interpreted as a definitive annotated or unannotated event. Annotation absence, unsupported chromosome, assembly mismatch and query failure remain distinct.

Annotation source/query metadata and tables are exported separately from the frozen RNA result. Display prioritization and transcript caps do not alter RNA counts or the annotation matching universe. No annotation-dependent redefinition of depth, junction support, genotype groups, PSI or causal inference is introduced.

## Local splice sequence view

The RNA allele evidence figure uses existing group `mean_depth`, `analyzed_n`, `count_sum` and `count_mean` without recomputation or denominator changes. All called genotypes use the same genomic coordinates and depth scale. Failed/unmeasured evidence remains unavailable; it is not drawn as zero. The existing RNA base observation table and its CIGAR/base-quality exclusions are retained.

A plus-strand donor motif comprises the first two intronic bases, and its acceptor motif the last two. On the negative strand these roles reverse and the genomic sequence is reverse complemented before comparison with GT/AG. Whole transcript exon models determine annotation boundaries; observed RNA junction boundaries supply the measurement evidence. Motif disruption or creation means the selected DNA substitution changes a two-base canonical motif under these assumptions. It does not establish transcript expression, phased allele-specific splicing or causality, and does not reconstruct other variants on a haplotype.

Reference-mismatched and non-SNV alleles do not receive an inferred sequence effect. Canonical creation outside an observed/annotated boundary is sequence-only evidence and is explicitly distinguished from junction support. Sequence source identity, genomic coordinates, transcript strand, REF/ALT and the display interval accompany the splice evidence exports.
