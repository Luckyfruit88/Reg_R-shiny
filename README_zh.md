# RegTools Shiny MVP

R/Shiny 网页参数界面 + 原生 RegTools junction 提取 + samtools depth + R 自定义计数。

**交付状态：代码和测试用例已编写；交付环境没有 R / samtools / RegTools，因此尚未执行 R 语法解析、R 单元测试、原生工具集成测试或浏览器测试，也没有对你的真实 BAM 验证。附带的 Python oracle 只独立核对模拟数据的预期值，不验证 R 程序本身。请先跑模拟测试，再用于真实样本。**

## 1. 适用范围

按经过 splice-aware 基因组比对的常规短读长 RNA-seq BAM 设计。需要 BAM 及匹配的 BAI/CSI 索引；坐标须与 BAM 的参考组装、染色体名称一致。基因组 DNA BAM 或转录本坐标 BAM 不能直接当作这种输入。

推荐本机或可信 Linux 内网服务器。RegTools 官方支持 Linux/macOS；Windows 可将计算端置于 Linux 环境。此原型不含公共服务所需的认证和多租户隔离，默认只监听 127.0.0.1。

BAM 只读；中间数据写入独立临时目录，正常退出或出错后清理。不自动改写、排序或重新索引源 BAM。进程被强制终止时可能残留临时文件，管理员需按本地数据保留政策清理。

## 2. 安装与启动

在已经安装 conda 的 Linux/macOS 终端，进入解压后的项目目录：

```bash
conda env create -f environment.yml
conda activate regtools-shiny
Rscript tests/run_tests.R
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838, launch.browser=TRUE)'
```

项目内 `environment.yml` 是待解析的依赖清单，不是已经在本交付环境成功安装的锁文件。若发生依赖求解失败，先检查 conda-forge / bioconda 的访问与平台支持。安装成功后可保存实际版本：

```bash
conda list --explicit > environment.explicit.txt
```

初次网页使用保留“内置模拟数据”，点击“读取区间 / 应用过滤参数”。无需提供任何个人测序数据。

### 接入你的 BAM

退出当前应用，然后在同一个 conda 环境设置 BAM 所在目录并重启：

```bash
export REGTOOLS_BAM_DIR="/absolute/path/to/bam_directory"
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838, launch.browser=TRUE)'
```

目录在启动时扫描一次，不递归；每个 BAM 及其 `.bam.bai`、`.bai`、`.bam.csi` 或 `.csi` 索引放在同一目录。解析后指向目录外的符号链接不会列入白名单。没有设置环境变量时使用项目的 `data/`。

只在源数据尚未正确排序和索引时，手动生成新的输出文件，且先确认输出文件名不存在：

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

网页上选择“服务器目录”和 BAM，输入染色体及固定统计区间，点击读取。再进入“Junction 表”，点击一行作为目标；也可自行输入目标内含子坐标。选择真实 BAM 后不要继续使用 `chrDemo`。

远程服务器仍可只绑定 127.0.0.1，通过你已有的 SSH 转发策略访问；不要为方便访问而直接将未经认证的原型暴露到公网。

## 3. 五项指标的明确口径

这些是本原型的应用定义，不声称它们是 RegTools 的五个原生输出列。

| 指标 | 原型中的含义 |
|---|---|
| `junction_reads` | 固定区间内，支持至少一个通过 RegTools 的 junction 的 alignment 记录数；对 read_id 取并集。不是所有 junction score 相加。 |
| `read depth` | samtools depth 的逐碱基结果；主卡片显示固定统计区间的零覆盖在内的平均值，同时导出全部位置。 |
| `junction_exact_count` | 内含子左右边界均与目标完全一致的支持记录数；合并所有链。 |
| `junction_near_count` | 左右边界分别都在目标 ±δ bp 内、但不是 exact junction 的支持记录数；对 read_id 取并集。 |
| `junction_per_100_reads` | `100 * junction_exact_count / denominator_reads`。分母为 0 时为 NA，不伪装成 0。 |

### 分母

`denominator_reads` 是经过相同 alignment 过滤后，至少一个 M、= 或 X 比对块真正覆盖固定分析区间的记录数。仅 reference span 与该区间重叠、实际只是 N 或 D 跨过区间的记录不计入。

按 **primary alignment 记录** 计数，paired-end 两端分别计数；不是 fragment、QNAME 去重计数或 UMI 数。记录层面的并集只防止一个多 junction read 在同一汇总项中重复计数，不进行分子去重。

精确目标支持数与其 RegTools score 按同一链口径对齐。对于多个 near junction，先取 read 集合再计数，避免简单加总重复计数。一个多 junction read 理论上可能同时支持 exact junction 和另外一个 near junction，因此两项未必互斥；`exact_or_near_reads` 专门导出二者并集，不应盲目相加。

本定义不是 PSI，也不是“每 100 个覆盖该单一剪接位点的 reads”，不是每 100 个全库 reads。更换分母定义会更换指标含义。多样本比较须固定区间、过滤条件、计数单位、坐标体系和链口径。

### Near

这里默认 **两端同时接近**，δ 默认 5 bp；不是只要 donor 或 acceptor 任意一端接近，也不是变异位点周围任意 junction。适合检查边界轻微偏移，但不能完整代表共享 donor/acceptor 的所有替代剪接。

界面 δ 为 0–50 bp；后端可接受 0–1000 bp。目标及 ±δ 搜索范围必须完全处于已读取区间，并且两侧至少各留 1 bp，否则报错要求扩大读取区间，不静默返回被截断的 near count。

### Depth

Depth 以参考位置逐碱基计算，M、=、X 参与；N 跳过，D 默认不计覆盖。未开启 samtools depth 的 `-s`，因此 paired-end 重叠部分按两个 reads 计，与 alignment 计数单位保持一致。

BaseQ 仅影响 depth，并不自动成为 RegTools 的逐 read anchor 碱基质量阈值。主卡片平均值包括整个统计区间的内含子和零覆盖位置，不应叫作平均外显子深度；还导出目标左右各最多 50 bp 窗口的均值，边界处会裁剪到统计区间。这些窗口未经过 GTF 验证，不能保证全是外显子。

## 4. 坐标约定

界面：1-based、两端均包含的**被跳过内含子**，不是两个外显子的边界碱基。

例如 `POS=181, CIGAR=20M100N30M`：

```text
左侧比对块       181–200
内含子           201–300  ← 界面输入
右侧比对块       301–330
内部 BED 坐标    [200, 300)
```

RegTools BED12 的 `chromStart/chromEnd` 包含两侧 anchor，不能直接作为精确内含子坐标：

```text
intron_start0 = chromStart + 第一个 blockSize
intron_end0   = chromEnd   - 第二个 blockSize
intron_start1 = intron_start0 + 1
intron_end1   = intron_end0
```

负链也保存染色体数值由小到大的左右端；不要直接把较小坐标称为 donor。原型保留各 junction 的 strand 列，但五项汇总使用所有链（包括 `?`），尚未实现链特异的分母。

## 5. 过滤与原生计数一致性

初始参数 MAPQ≥20、BaseQ≥0、anchor≥8、intron 70–500000 bp。这只是原型起点，尤其 MAPQ、最短内含子、文库方向需要按比对器、物种和实验设计调整。

固定排除：UNMAP(4)、SECONDARY(256)、QCFAIL(512)、SUPPLEMENTARY(2048)，合计 `-F 2820`。已标记 duplicate 默认保留，可显式排除（mask 变为 3844）；代码不进行 duplicate 标记，也不推断哪些重复应被去除。

MAPQ 并非跨比对器一致的“唯一比对”标签；255 表示 mapping quality 不可用。`NH==1` 可选，但缺失 NH 时也会被排除，界面会提醒。未提供 UMI 模式。

RegTools 的 anchor 规则是 junction 层面的：左右最小 anchor 条件可以由不同 reads 满足。代码不再擅自加一个“每条 read 两端都 ≥8 bp”的过滤，否则会改变原生 count。RegTools 会依据自身 CIGAR 规则处理 anchor；不能把通用 M CIGAR 的 anchor 简化宣称为经过独立错配验证的完美匹配。

每次计算后将 accepted junction 的逐 read CIGAR 证据计数，与相同输入、相同链模式下 RegTools BED12 score 逐项核验。不同则中止，不输出看似成功的统计。这只是实现一致性检查，不证明比对或生物学解释正确。

RF/FR 模式按 RegTools flags 逻辑匹配；XS 模式缺失标签保留 `?`。原型不使用 FASTA intron-motif 推断。需要已知方向对照验证自己的文库，单端数据尤其不能直接套用 paired-end 文库假设。

## 6. 交互与资源限制

重操作只在“读取区间”点击时进行；Shiny ExtendedTask + future worker 执行一次参数快照。目标、δ、显示窗口变化复用当前已完成数据，不自动再扫描 BAM。更改输入 BAM、统计区间、MAPQ、anchor 等后，界面显示旧结果提示，重新读取后才应用。

只保留本会话最近一次分析结果，不是按 BAM/参数哈希持久化的跨任务缓存。图形缩放不改变统计分母。默认 2 个 future workers；生产环境需按服务器核数、内存与并发用户重新配置。

上限：单次分析 250 kb、100,000 条候选 alignment、局部 SAM 256 MiB。单个外部命令超时为 600 秒。这些是原型保护阈值，不是运行耗时承诺；超过阈值会明确报错，不抽样。未实现任务取消、跨会话作业队列或完整磁盘配额管理。

图形中 junction 弧线只显示当前视窗完全包含且支持数最高的 30 条；计数和导出不受这个展示上限影响。弧线不是完整 Sashimi/transcript annotation。

## 7. 输出

- `junction_metrics.csv`：五项指标、分母、坐标、两侧窗口 depth 和并集计数。
- `junctions.csv`：native BED12 各列、内含子 0/1-based 坐标、score、链和目标分类。
- `depth.tsv`：固定区间所有碱基的位置与覆盖度，包括零覆盖。
- `junction_read_evidence.csv`：所有 accepted junction 的逐 read 证据，包括 QNAME、FLAG、MAPQ、CIGAR、NH、XS；不导出 SEQ/QUAL。
- `parameters_and_audit.json`：本次参数、目标定义、BAM 大小/mtime、工具版本、命令日志、核验状态。文件可能含样本名、QNAME 或本地路径，分享前请检查隐私。

文件大小/mtime 是轻量身份检查，不是 BAM 的密码学校验。中间 BAM/日志路径在计算结束后失效；JSON 留存的是运行记录，不是可不加修改重放的 shell 脚本。正式可复现分析还应固定代码、依赖版本与输入文件校验值。

## 8. 模拟测试与验收标准

默认 chrDemo:101–450，目标 intron=201–300，δ=5、MAPQ≥20：

| 项目 | 预期 |
|---|---:|
| denominator_reads | 4 |
| junction_reads | 3 |
| junction_exact_count | 2 |
| junction_near_count | 1 |
| junction_per_100_reads | 50 |
| mean_read_depth | 200/350 ≈ 0.5714286 |

模拟数据同时包含低 MAPQ、secondary、supplementary、只以 N 跨过区间的反例记录。把 δ 从 5 改为 1，near 应变为 0；改为 2 则为 1。把 MAPQ 改为 61 后重新读取，分母为 0、rate 为 NA。图形缩放不能改变任何统计值。

`tests/run_tests.R` 含 CIGAR、坐标、文库方向、模拟数据、原生工具一致性、零分母和 near 边界测试。**该 R 测试脚本在交付环境未执行**。`tests/oracle.py` 是无第三方依赖的独立期望值核验，在交付环境已执行；不应拿它替代 R 集成测试。

在真实 BAM 上，先选择一个熟悉的小区域，核对精确坐标、链模式及同一过滤条件下的命令行 RegTools 输出；再开展批量使用。不要直接拿另一工具的默认显示数比较，却忽略 duplicate、MAPQ、paired-end、anchor 等口径差异。

## 9. 暂未实现

多样本队列、GTF/FASTA 注释、已知/新 junction 分类、fragment/UMI 计数、CRAM、单细胞分组、链特异分母、正式 PSI、共享 donor/acceptor 事件分析、原始 supporting BAM/FASTQ 下载、生产级认证权限与任务取消。

## 10. 一手资料（2026-09-16 查阅）

RegTools 命令、BED12 score/坐标与 anchor 定义：
https://regtools.readthedocs.io/en/latest/commands/junctions-extract/

RegTools 1.0.0 源码（链推断、CIGAR、计数）：
https://raw.githubusercontent.com/griffithlab/regtools/1.0.0/src/junctions/junctions_extractor.cc

RegTools 官方仓库 / 系统要求：
https://github.com/griffithlab/regtools

Bioconda 安装与版本配方：
https://bioconda.github.io/recipes/regtools/README.html

samtools view：
https://www.htslib.org/doc/samtools-view.html

samtools depth：
https://www.htslib.org/doc/samtools-depth.html

SAM 格式与 flags：
https://www.htslib.org/doc/sam.html
https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex

UCSC BED 坐标规范：
https://genome.ucsc.edu/FAQ/FAQformat.html#format1

Shiny 非阻塞任务：
https://shiny.posit.co/r/articles/improve/nonblocking/

processx 的参数向量调用和超时：
https://processx.r-lib.org/reference/run.html
