# Reg_R-shiny

R Shiny 版 RegTools 局部 RNA-seq junction 分析原型。

通过网页设置 BAM、固定分析区间、过滤条件、目标 junction 和 near 容差，查看 junction reads、read depth、junction exact count、junction near count 与 junction per 100 reads，并导出结果和逐 read 证据。

## 使用说明

完整安装步骤、统计口径、坐标约定、测试和限制见 [中文说明](README_zh.md)。

在已安装 conda 的 Linux/macOS 环境中，从项目根目录运行：

```bash
conda env create -f environment.yml
conda activate regtools-shiny
Rscript tests/run_tests.R
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838, launch.browser=TRUE)'
```

默认使用内置模拟数据。接入真实 BAM 前，将 `REGTOOLS_BAM_DIR` 设置为 BAM 及其匹配索引所在目录，再启动应用。不要将真实测序数据或个人路径配置提交到仓库。

## 验证状态

- 原型包文件 SHA-256 校验通过；`SHA256SUMS.txt` 保留原交付文件清单。
- 独立 Python 模拟数值核验通过：`python tests/oracle.py`。
- R 语法解析、R 单元测试、samtools/RegTools 集成测试及浏览器测试尚未执行。
- 未使用真实 BAM 验证；此原型不是经过临床验证的软件。

源码上传不等于应用部署。该原型默认仅监听 `127.0.0.1`，没有公共服务所需的登录认证或多租户隔离。

## 项目文件

- `app.R`：Shiny 界面与交互。
- `backend.R`：BAM 读取、RegTools 提取、depth 计算和指标汇总。
- `environment.yml`：待安装验证的 conda 依赖清单。
- `tests/`：R 测试、Python 期望值核验和验证状态。
- `README_zh.md`：完整中文文档及一手资料。
