# Copy Number–Driven Dependency Analysis

> 初始准备阶段的历史 README。当前模块化 R 框架、正式案例和数据状态见 [项目 README](../README.md)。

用于未来任意 Gene A / Gene B 的拷贝数、表达、CRISPR dependency、突变、癌种及基因组背景分析。本轮只准备数据库、QC 和数据访问接口，不执行生物学分析，也不创建完整分析 Skill。

## 本次准备状态

官方免验证码目录自动识别出的最新 DepMap Public release 为 **26Q1**。用户明确要求暂时跳过 DepMap 下载，因此该版本只是已识别、已选择；没有安装 DepMap 实体数据库。不能把版本检查的 “up to date” 理解为文件已下载。

TCGA 的实际下载、预处理、空间占用、QC 和所有未下载文件见 [准备报告](docs/DATABASE_PREPARATION_REPORT.md)。逐文件 URL、下载日期、字节数和 SHA256 见 [manifest](data/manifests/data_manifest.csv)，来源定义见 [DATA_SOURCES](docs/DATA_SOURCES.md)。

## 目录与版本

```text
data/
  raw/depmap/<release>/       原始文件；禁止编辑、覆盖或删除
  raw/tcga/pancan/           Xena 原始 gzip 文件
  processed/depmap/<release>/ModelID × genes，ZSTD Parquet
  processed/tcga/            按基因存储的 ZSTD Parquet + DuckDB views
  manifests/                下载清单、源数据字典、QC、验证结果
scripts/download/           下载和版本检查
scripts/preprocess/         分批预处理
scripts/utils/              校验和报告
scripts/data_access.py      统一查询接口
config/                     release 选择、文件白名单、TCGA 癌种映射
results/                    预留分析输出目录（本轮无分析结果）
logs/                       日志
docs/                       来源说明、完整报告
tests/                      基础设施安全与映射测试
```

DepMap processed 也按 release 分目录，避免更新时覆盖旧版本。`config/current_release.json` 选择当前版本，接口只读取该版本。没有 WGS CN 的 Model 可补用同 release 的 `OmicsCNGene.csv`；26Q1 官方目录不提供此文件，因此本次不从其他 release 补取。所有矩阵保留全部基因。

## Windows / Git Bash 运行

使用 C 盘已有的 Python 3.12；不要调用 `python` / `python3` Windows Store stub。依赖安装在项目 `.runtime/`，脚本自动加入此路径。下面的 `project_root` 是 Windows-native 路径，供 Windows Python 识别。

```bash
project_root='D:/CodexProjects/Copy Number–Driven Dependency Analysis'
cd '/d/CodexProjects/Copy Number–Driven Dependency Analysis'
/c/Python312/python.exe -m pip install --target "$project_root/.runtime" -r "$project_root/requirements.txt"

# 只检查 release；不会下载或替换数据库
/c/Python312/python.exe "$project_root/scripts/download/check_updates.py"

# 本次已按用户要求跳过 DepMap。将来明确决定下载时再使用此命令。
# 先在官方浏览器通过验证，导出具有 release,filename,url,md5_hash 的目录 CSV。
/c/Python312/python.exe "$project_root/scripts/download/download_depmap.py" --url-catalog 'D:/path/to/browser_verified_catalog.csv'

# 新 release 必须显式指定；旧 raw 和 processed 均保留
/c/Python312/python.exe "$project_root/scripts/download/download_depmap.py" --release '<new release>' --url-catalog 'D:/path/to/browser_verified_catalog.csv'

# 重新运行已验证文件会先检查 SHA256，然后跳过；不会覆盖 raw
/c/Python312/python.exe "$project_root/scripts/download/download_tcga.py"

# 仅下载某个 TCGA 文件
/c/Python312/python.exe "$project_root/scripts/download/download_tcga.py" --only 'mc3.v0.2.8.PUBLIC.xena.tsv.gz'

# 重新生成 processed
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_depmap.py"
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_tcga.py"
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_tcga_mutation.py"

# 完整文件校验和单基因读取与 raw 一致性检查
/c/Python312/python.exe "$project_root/scripts/utils/validate_database.py"
/c/Python312/python.exe "$project_root/tests/test_infrastructure.py"
/c/Python312/python.exe "$project_root/scripts/utils/build_report.py"
```

同一时间只运行一个下载命令，避免下载 manifest 的并发写入。下载中间文件使用 `.part`；断点续传前比较 URL、ETag 和长度。服务器不支持安全续传、文件远端身份变化或 checksum 不匹配时保留中间文件并报错，不自动清理成功文件。已有 raw 校验失败时拒绝覆盖，需要单独检查。

当前公开 MC3 MAF 的 GDC 历史链接失效，重跑完整 TCGA 下载命令会继续记录这个可选文件失败。Xena MC3 表是稳定的补充来源，但它保留原 MAF 的 `FILTER=PASS` 记录，不能当成完整原 MAF。

## 数据访问

分析代码统一导入 `scripts.data_access`，不用硬编码 CSV 路径：

```python
from scripts.data_access import (
    get_depmap_cn, get_depmap_expression, get_depmap_dependency,
    get_depmap_model_metadata, get_depmap_mutation,
    get_tcga_cn, get_tcga_gistic, get_tcga_expression, get_tcga_metadata,
    request_gene,
)
# 任意合法 gene symbol；这里只展示接口，不执行分析
cn = request_gene("TP53", database="tcga", dataset="cn")
meta = get_tcga_metadata()
# 概率矩阵使用 get_depmap_dependency(gene, probability=True)
# mutation 默认 damaging；kind="hotspot" 可读取热点二分类矩阵
```

DepMap 查询通过 Parquet column projection 只读 ModelID 和目标基因列。TCGA 每个基因单独占一个 Parquet row group，Python 接口通过独立基因索引直接定位并只读取一个 row group；缓存仅含基因索引和 Parquet footer，不缓存整个矩阵。`Values` 按独立 samples 表的顺序返回。DuckDB 文件仅保存可选 SQL views，不复制完整矩阵。移动项目目录后重跑 `preprocess_tcga.py` 可以重建 views；Python 接口本身使用当前项目路径。

所有数值保留源数据尺度与 float64 精度。TCGA Toil `tcga_RSEM_Hugo_norm_count` 是 **log2(norm_count+1)**，不是 TPM。CNA 的 `-2/-1/0/1/2` 分类仅来自 TCGA thresholded GISTIC，不将其直接套到 DepMap relative CN。患者 ID 从 TCGA barcode 前 12 字符提取，肿瘤/正常来自官方 sample-type code，癌种来自 phenotype 的明确疾病映射。多样本患者不会被自动合并。

DepMap 实体文件未下载时接口抛出明确的 `DataUnavailableError`。未知基因、同名多列和无法唯一映射 ModelID 的情况会报错，不自动猜测、平均或混用版本。

## 更新与可复现性

定期运行 DepMap release 检查；出现新版只提示，明确选择后才能下载。每个 release 独立保留 raw、processed 和清单。TCGA 视作固定历史 cohort，记录下载时间、源版本、URL、ETag、SHA256，不定期重新下载。源文件变化时保留旧数据另行决定更新。

`.gitignore` 排除了 `data/raw/`、`data/processed/`、本地依赖和日志；只提交脚本、配置、manifest 和文档。不要把大型数据库加入 Git。
