# Reg_Shiny：SCC 成员各自 clone、各自运行

Reg_Shiny 用 DNA VCF 的 genotype 分组，对照匹配 RNA-seq BAM 的 depth、junction 和局部序列证据。**默认连接 FHS WGS/RNA 数据，也可以在界面选择其他 SCC 文件。** 所有计算与测序数据均留在 SCC。

## 在 Interactive RStudio Server 启动

1. 在 SCC OnDemand 新建 **RStudio Server**，选择 **R 4.5.2、3 核、无 GPU** 和你有权限使用的项目，Extra qsub options 填 `-l mem_per_core=4G`。使用自己可写的持久化项目目录；首次 GENCODE 索引约需 1 GB，不要把缓存和作业结果放在节点临时目录。
2. 在 RStudio **Terminal** 执行：

   ```bash
   git clone https://github.com/Luckyfruit88/Reg_R-shiny.git
   ```

   如 GitHub 要求认证，请使用自己已有权限的账号。GitHub 代码访问和 SCC 数据权限是两回事。
3. 通过 RStudio **File → Open Project** 打开 clone 中的 `Reg_R-shiny.Rproj`。仅在 Terminal 中 `cd` 不会改变 R Console 的工作目录。
4. 在 **R Console** 执行：

   ```r
   source("run_app.R")
   run_reg_shiny()
   ```

   也可以打开 `app.R` 后点击 **Run App**。启动器会把缺失的 SCC 原生工具路径载入当前 R 会话，使用 `samtools/1.23 regtools/1.0.0 bcftools/1.23`。R 4.5.2 的共享库已提供所需 R 包；其他环境缺包时可执行 `source("scripts/install_dependencies.R")`，仅安装到个人 R library。

首次启动会在 **Data sources** 自动准备 FHS：从原始对应表建立本用户的 VCF/BAM 清单，并建立 GENCODE v48 区域索引；界面显示当前阶段。准备完成后进入 WGS 页面，后续复用经过身份校验的本地缓存。无需复制 Jinjie 的配置、索引、作业目录或会话链接，也无需额外设置 `R_ENVIRON_USER`。

## 选择自己的 BAM / VCF

在 **Data sources → Custom SCC files** 中：

- 填写数据集名称及真实参考组装（如 GRCh38）。选择一个带索引的 VCF/BCF，或填写按 chromosome 映射到 VCF 的 TSV。
- 指定 BAM 文件/目录。进行 genotype 对照时，提供 `vcf_sample` 与 BAM 绝对路径的对应表；或选择严格匹配 BAM read-group `SM` 与 VCF sample ID。程序不猜测文件名与受试者的关系。
- 可填写匹配的参考 FASTA 和 GENCODE v48 GTF/已有 SQLite 索引。GENCODE v48 与当前 GT/AG 序列叠加需要兼容的 GRCh38；其他组装不会强行套用该注释。
- 点击 **Prepare and connect**。成功后才切换数据源，旧数据的查询和图不会被当成新数据；失败时保留原连接。

只有 BAM 时可以不填 VCF，使用 **Single BAM / demo**。[对应表与 VCF registry 示例](docs/DATA_SOURCES.md) 均为虚构样例，不含真实个体数据。

## 查询和查看结果

输入 `CHROM:POS` 或 `CHROM:POS:REF:ALT`，选择对应记录后运行 **All matched BAMs**。界面报告实际配置、选中、成功和失败样本数（FHS 对应唯一参与者）；不存在的 genotype 不会补造，缺失 genotype 不是 0/0。

批量查询时，在 **Batch variants** 中每行输入一个坐标，然后点击 **Submit independent 16-core jobs**。当前窗口大小和过滤条件分别围绕每个 variant 冻结。查询与提交在后台处理，已提交的任务可在后续条目准备时开始计算。每条输入都有独立的提交结果或失败原因；同一坐标有多个 REF/ALT 时，需要改用明确的 `CHROM:POS:REF:ALT`。在任务列表中打开某个任务，即可查看它自己的结果和日志。

**RNA allele evidence** 保留 RNA 碱基计数表，并按 genotype 显示逐碱基平均覆盖度、参考序列、GT/AG 改变和实测 junction 支持。不同 genotype 共用刻度；测量为 0 与不可用 NA 分开。GENCODE 页显示参考 transcript/exon/CDS/UTR 和精确 junction 注释。上述结果用于检查剪接证据，不单独证明因果，也不是正式 PSI 或关联检验。

全量比较作为独立 SCC 作业运行，关闭浏览器不会停止。回到相同的 prepared dataset，可在保存的作业列表中打开或恢复。原始 BAM/VCF 只读，不自动重建源索引。

## 每位成员的文件与权限

代码目录下 `.reg_shiny/uid-<uid>/` 保存本用户的私有配置、清单、GENCODE 索引与作业结果，均被 Git 排除。每个数据集有自己的作业存储，同一数据集可同时运行多个 variant 任务。每个任务独立保存准确的 VCF 记录、genotype/BAM 匹配、参数、代码快照、检查点、进度、日志与结果；某个任务失败不会停止其他任务。

默认 FHS 源数据需要已有 `mtdna-alcohol` 和 `sequencing` 权限。clone 代码不会新增数据权限。不要使用他人的私有 `.git`、缓存、用户作业目录或 OnDemand rnode 链接。

每个新全量任务独立申请 **16 核**（15 个 RNA worker 和 1 个协调进程）、4 GB/核，默认 12 小时上限。例如，3 个 variant 同时运行会占用 **48 核**，不会共享一份 16 核。实际并行启动时间取决于 SCC 排队和资源可用性；增加核数不保证耗时减半。

计费项目与时限可在启动前设置 `REGSHINY_SGE_PROJECT`、`REGSHINY_FULL_WALLTIME`。旧的 `REGSHINY_FULL_CORES` 不再改变新任务的 16 核配置；历史任务恢复时保留其冻结的原资源配置。交互界面仍按 RStudio 核数使用 1–2 个后台 worker，不承担全量计算。

更新时先停止 app，在本 clone 的 Terminal 中执行 `git pull --ff-only`，然后重新启动。更新不会主动删除私有缓存和作业；重连时仍检查源文件身份。

[测试与验收说明](docs/VALIDATION.md) · [完整指标定义](docs/METRICS.md)

## 许可证

Reg_Shiny 的源代码和文档采用 [MIT 许可证](LICENSE)。上游来源与许可声明见 [第三方声明](THIRD_PARTY_NOTICES.txt)。

软件许可证不授予 FHS 或其他受控数据、BAM/VCF 文件及外部参考资源的访问或再分发权；这些资源仍遵循各自的访问与使用条款。第三方软件依赖保留各自的许可证。
