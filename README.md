# Reg_Shiny

Reg_Shiny compares RNA-seq splice evidence by DNA genotype at a selected WGS variant. Run your own clone inside **SCC Interactive RStudio Server**. It connects to the FHS WGS/RNA sources by default, and the **Data sources** tab lets you select other indexed VCF/BCF and BAM files already on SCC.

[中文使用说明](README_zh.md) · [Data-source formats](docs/DATA_SOURCES.md) · [Metric definitions](docs/METRICS.md)

## Start in SCC RStudio

1. Start an SCC OnDemand **RStudio Server** session with **R 4.5.2**, **3 cores**, no GPU and your authorized project. Set Extra qsub options to `-l mem_per_core=4G`. Choose your own writable, persistent SCC project directory; the first GENCODE index needs approximately 1 GB. Do not use another member's clone or a node's temporary directory for persistent results.
2. In the RStudio **Terminal**, clone the repository using your existing GitHub access:

   ```bash
   git clone https://github.com/Luckyfruit88/Reg_R-shiny.git
   ```

   If GitHub requests authentication, use your existing authorized account. Repository access and SCC data-group permissions are separate.
3. In RStudio, use **File → Open Project** to open the cloned `Reg_R-shiny.Rproj`. This sets the R Console's working directory correctly; a Terminal `cd` alone does not change it.
4. In the R Console, run:

   ```r
   source("run_app.R")
   run_reg_shiny()
   ```

   Opening `app.R` and using **Run App** uses the same initialization. The launcher loads missing SCC native-tool paths into the R session. Supported modules are `samtools/1.23 regtools/1.0.0 bcftools/1.23`; SCC R 4.5.2 already supplies the R packages. If R dependencies are missing in another environment, `source("scripts/install_dependencies.R")` installs only missing packages into your personal library.

On first launch, **Data sources** automatically prepares the default FHS connection and shows preparation progress. It builds private manifests from the original source maps and, when available, a private GENCODE v48 index. Reopening a prepared connection reuses its validated files. The app then opens WGS genotype comparison. No personal deployment path, copied participant manifest, manually set `R_ENVIRON_USER`, or another user's rnode link is required.

Default FHS access requires existing **mtdna-alcohol** and **sequencing** authorization. Cloning code does not grant data access. The application does not change source-file permissions, rewrite BAM/VCF files or create their indexes.

## Choose data and inspect a variant

- Keep **FHS on SCC (default)**, or select **Custom SCC files** in **Data sources**. Enter a single indexed VCF/BCF or a chromosome-to-VCF registry, then a BAM directory/file or an explicit DNA-to-RNA mapping TSV. For genotype comparison, use an explicit mapping or verified exact BAM read-group `SM` / VCF sample-name matches. Filenames are never guessed as participant IDs. See the [format examples](docs/DATA_SOURCES.md).
- Click **Prepare and connect**. Validation succeeds before the active dataset changes. A failed preparation preserves the current connection; a successful switch clears the previous viewer and keeps persistent jobs under their original dataset.
- Search `CHROM:POS` or `CHROM:POS:REF:ALT`, select the exact VCF record and compare genotypes. **All matched BAMs** is the default. Optional preview is labelled and limited only when explicitly selected. A missing genotype is not invented or treated as reference homozygous.
- Inspect depth, exact junction support, GENCODE transcript/exon/CDS/UTR models and the base-resolution **RNA allele evidence** view. The latter keeps the RNA base table and aligns genotype mean coverage with reference sequence and GT/AG changes. Raw counts and genotype differences are descriptive evidence, not formal PSI or causal inference.

VCF may be left blank for **Single BAM / demo** use. Optional FASTA/GENCODE inputs are configurable. GENCODE v48 and the current splice-motif reference overlay require compatible GRCh38 inputs; other declared assemblies can still use their compatible BAM/VCF counting views without a false GRCh38 overlay.

## Private state and full-cohort jobs

Each clone creates `.reg_shiny/uid-<uid>/` with private profile directories. These hold source manifests, source identities, annotation indexes and results and are excluded from Git. Original sequencing files stay in their original authorized locations. A profile is bound to its source identity; source changes require preparing the connection again.

Full comparisons submit a separate SCC job, stream per-sample contributions and save checkpoints. Closing the browser does not cancel that job. Reopen or resume it from the same prepared dataset's saved-job list. Each dataset's job store admits one active full comparison; separate profiles/clones are separate stores, so users must account for their total active jobs. Defaults are 8 cores, 4 GB/core and a 12-hour walltime. In the R Console, before starting the app, you can set `REGSHINY_SGE_PROJECT`, `REGSHINY_FULL_CORES` and `REGSHINY_FULL_WALLTIME` for your authorized allocation.

The interactive session uses two background workers when `NSLOTS >= 3`, otherwise one, while the main process handles the UI. `REGSHINY_UI_WORKERS` may reduce this within the detected budget. Thread limits set after RStudio has started cannot retroactively reconfigure already initialized numeric libraries.

To update, stop the app, run `git pull --ff-only` in the repository's Terminal, then restart. Do not commit `.reg_shiny`, participant tables, sequencing files or analysis outputs. Existing private data and job directories are not deleted by an update; frozen-source checks determine whether a profile/result can be reused.

## Validation and scope

The repository includes synthetic native-tool, backend, source-import, asynchronous UI and saved-job tests. [Validation instructions](docs/VALIDATION.md) distinguish those checks from real SCC acceptance. No controlled dataset or individual-level output is distributed with the code. An owner-account run and permission inspection do not substitute for a member's own GitHub/SCC login.

BAM counts use retained primary alignment records (paired ends separately), with explicit index, coordinate and reference checks. NA and measured zero remain different. RNA base observations do not replace DNA genotype calls; inferred GT/AG changes apply only to the selected SNV and do not reconstruct an allele-phased RNA haplotype. See [METRICS.md](docs/METRICS.md) for complete definitions.

GENCODE source: [human release 48](https://www.gencodegenes.org/human/release_48.html). Native tools: [RegTools](https://github.com/griffithlab/regtools), [samtools](https://www.htslib.org/), [bcftools](https://samtools.github.io/bcftools/).
