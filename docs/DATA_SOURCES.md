# Data sources

All paths refer to files on the host running Reg_Shiny (SCC for the documented workflow). Selecting a path does not upload, copy, modify or index the original BAM/VCF. Use the same assembly and literal chromosome names across the inputs; there is no automatic liftover or chromosome-prefix conversion.

## Default FHS connection

The FHS preset uses the original RNA-ID / BAM-location CSVs under the project's shared FHS RNA-seq `IDs` directory and the sequencing project's `freeze_10a/passgt.minDP10.fhsids` chromosome VCFs. It reconstructs each member's manifests from those original sources. The exact default paths are centralized in `default_fhs_spec()` in `data_sources.R`.

The preset includes chr1–22 and chrX and all unique DNA/RNA matches with a readable BAM and existing index. It does not use the old 8,973-pair candidate list or filter to a covariate-selected research cohort. Preparation records exclusions and actual counts privately. Every member still requires the source owners' existing SCC data permissions.

GENCODE v48 is indexed once in this account's private resources and reused across compatible profiles. No index or sample manifest in another user's Reg_Shiny deployment is required. Reference/annotation availability is reported independently of RNA-count availability.

## Custom WGS + RNA inputs

Choose one DNA input in **Data sources → Custom SCC files**:

- An indexed `.vcf.gz`, `.vcf.bgz` or `.bcf` containing one or more chromosomes. Automatic configuration uses contigs with records reported by the existing index, so unused header declarations do not block a normal chromosome file.
- A tab-separated registry with columns `chrom`, `vcf`, `build`, with one row per unique literal chromosome. One multichromosome VCF may appear in several rows. VCF paths must be absolute. A readable `.tbi` or `.csi` must already exist beside each DNA file. See [vcf_registry.example.tsv](../examples/vcf_registry.example.tsv).

Choose one matching rule:

1. **Explicit mapping TSV** (recommended): `vcf_sample` is an exact VCF sample-column name; `bam` is an absolute BAM path. Optional `build` must agree with the selected assembly. See [sample_mapping.example.tsv](../examples/sample_mapping.example.tsv). Participant/RNA library IDs that differ require this explicit mapping.
2. **Exact read-group SM matching**: specify a BAM file or directory. Each BAM must have one unambiguous read-group sample identity matching the VCF sample name exactly. The importer does not infer identity from a filename, prefix or numeric conversion.

Every mapped DNA sample and BAM must be unique. Duplicate physical BAMs, multiple BAMs for one DNA ID, ambiguous read-group identities and incompatible headers require correction. The app never resolves these cases by selecting an arbitrary row. Reads from technical replicates must be resolved by an explicit upstream data decision.

BAM directory discovery examines direct `.bam` children; it is not a recursive filesystem scan. BAM indexes may use `.bam.bai`, `.bai`, `.bam.csi` or `.csi`. Other sample rows in a VCF can remain unmatched, and the interface reports the resulting unmatched-DNA counts. Missing/partial DNA GT stays missing and does not join a reference-homozygous group.

## BAM-only and optional references

For Single BAM exploration, omit both DNA fields and the sample mapping, and specify one BAM or a directory. A synthetic `chrDemo` remains available independently of source preparation.

An optional FASTA needs an existing `.fai`. An optional GENCODE input is either a complete v48 GTF/GTF.gz to index or a compatible SQLite index produced by `scripts/build_annotation.py`. GENCODE v48 and the current splice-sequence overlay require GRCh38; other explicitly declared assemblies retain compatible DNA/RNA count views and do not receive a false GRCh38 annotation.

## Reuse and switch

**Prepare and connect** publishes validated manifests atomically. A failed attempt leaves the current connected dataset unchanged. A successful switch uses a fresh viewer, and results from the prior source are not reused under new labels. Each profile's saved jobs retain their own data identity.

Private state is stored at `.reg_shiny/uid-<uid>/` in this clone. Do not copy or edit another person's profile directories. Source file identity and prepared manifest hashes are checked when reconnecting. Changed inputs require preparing a new profile; the old job history is preserved.

To choose a custom initial profile in the R Console, copy [profile.example.json](../config/profile.example.json) to a private location, edit its server paths, and run:

```r
Sys.setenv(REGSHINY_PROFILE_CONFIG = "/absolute/path/to/my-profile.json")
source("run_app.R")
run_reg_shiny()
```

No personal manifest or reference index should be committed to Git. Only explicitly synthetic templates are tracked in `examples/`.
