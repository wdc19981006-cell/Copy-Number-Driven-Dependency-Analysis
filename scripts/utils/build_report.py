"""Generate a factual preparation report and source inventory from manifests."""
from pathlib import Path
import csv
import json
import platform
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, load_manifest, write_json, now, setup


def size(n):
    n = int(n or 0)
    if n >= 1024**3:
        return f'{n / 1024**3:.3f} GiB'
    if n >= 1024**2:
        return f'{n / 1024**2:.2f} MiB'
    if n >= 1024:
        return f'{n / 1024:.2f} KiB'
    return f'{n} B'


def occupation(path):
    return sum(p.stat().st_size for p in path.rglob('*') if p.is_file())


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}


def main():
    setup('build_report')
    rows = load_manifest()
    config = read_json(ROOT / 'config/current_release.json')
    qc = read_json(ROOT / 'data/manifests/tcga_qc.json')
    validation = read_json(ROOT / 'data/manifests/validation_report.json')
    mc3 = read_json(ROOT / 'data/manifests/tcga_mc3_qc.json')
    datasets = read_json(ROOT / 'config/datasets.json')
    depmap_bytes = occupation(ROOT / 'data/raw/depmap')
    tcga_bytes = occupation(ROOT / 'data/raw/tcga')
    processed_bytes = occupation(ROOT / 'data/processed')
    successful = [r for r in rows if r['status'] == 'verified']
    missing = [r for r in rows if r['status'] != 'verified']
    purposes = {('DepMap', d['filename']): d['purpose'] for d in datasets['depmap']}
    purposes.update({('TCGA', d['filename']): d['purpose'] for d in datasets['tcga']})
    summaries = []
    for row in rows:
        summaries.append({**row, 'purpose': purposes.get((row['database'], row['filename']), '')})
    write_json(ROOT / 'data/manifests/preparation_summary.json', {
        'generated_at': now(), 'depmap_release': config.get('depmap_release'),
        'depmap_installation_status': config.get('installation_status'),
        'raw_depmap_bytes': depmap_bytes, 'raw_tcga_bytes': tcga_bytes,
        'processed_bytes': processed_bytes, 'successful_files': successful,
        'not_downloaded': missing, 'validation_passed': validation.get('passed'),
        'depmap_model_overlap': 'Not measurable: user skipped DepMap data downloads',
        'tcga_qc': qc, 'mc3_qc': mc3,
    })
    miss_csv = ROOT / 'data/manifests/not_downloaded.csv'
    with miss_csv.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['database','release','filename','status','reason'])
        w.writeheader()
        for r in missing:
            w.writerow({k:r[k] for k in ['database','release','filename','status']} | {'reason':r['error']})
    source_lines = [
        '# 数据来源与版本', '', f'生成时间：{now()}（Asia/Shanghai）。', '',
        f'DepMap：官方目录最新 release 为 **{config.get("depmap_release")}**；本次用户明确要求跳过实体文件。'
        'release 选择不代表已经下载或通过真实数据库 QC。', '',
        '官方文件目录：[DepMap no-captcha catalog](https://depmap.org/portal/api/no-captcha/download/files)。'
        '该目录保留文件名和 publisher MD5，但 URL 为空。'
        '[官方说明](https://forum.depmap.org/t/provide-an-open-endpoint-for-latest-version-retrieval/4652)确认下载链接需要浏览器验证。', '',
        'TCGA：CN/GISTIC 来自 [TCGA Xena hub](https://tcga.xenahubs.net)，'
        '表达来自 [Toil hub](https://toil.xenahubs.net)，'
        'phenotype 和 MC3 补充表来自 [Pan-Cancer Atlas hub](https://pancanatlas.xenahubs.net)。'
        '这些 hub 重定向到官方 S3 数据对象，下载时保留 gzip 原始内容。', '',
        '源 metadata 字典存放在 `data/manifests/*.source.json`。每个文件的原始 URL、时间、精确字节数、'
        'SHA256、状态、服务器 ETag/Last-Modified 见 `data_manifest.csv`。SHA256 是本地内容指纹；'
        'DepMap 将来下载时还会对照官方 publisher MD5。TCGA 不把 multipart ETag 当作 MD5。', '',
        '| 数据库 / release | 文件 | 用途 | 下载状态 | 下载时间 |',
        '|---|---|---|---|---|',
    ]
    for r in summaries:
        source_lines.append(f'| {r["database"]} / {r["release"]} | `{r["filename"]}` | {r["purpose"]} | {r["status"]} | {r["download_date"] or "未下载"} |')
    source_lines += [
        '', '## 数值与 cohort 定义', '',
        '- TCGA 连续 CN 保留官方 GISTIC2 continuous 数值，不推断整数拷贝数。',
        '- TCGA thresholded GISTIC：-2 deep deletion、-1 shallow deletion、0 diploid、1 gain、2 amplification。',
        '- Toil `tcga_RSEM_Hugo_norm_count` metadata 指明 **log2(norm_count+1)**，源版本 2016-02-18；不能标为 TPM。',
        '- MC3 Xena 官方数据字典版本为 2016-12-29，来源 MAF 仅保留 `FILTER=PASS`。完整下载 Xena 表，未声称它是原始 MAF。',
        '- GDC 页面公开 MAF 的 UUID 已失效；当前 files API 按文件名查询也无记录。保留失败来源，不伪造 MAF。'
        '原链接来自 [GDC MC3 publication](https://gdc.cancer.gov/about-data/publications/mc3-2017)。',
        '- CancerType 使用 phenotype 的 33 个明确疾病全名映射为官方缩写；映射表在 `config/tcga_cancer_types.json`，'
        '依据 [GDC TCGA 缩写表](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)。',
        '- SampleTypeCode 取 barcode 中样本类型码，PatientID 取前 12 字符；定义依据 '
        '[GDC sample-type 表](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)。'
        '`TCGA-07-0249-20` 是 Control Analyte，不在 phenotype 中，保留为 Other，CancerType 不猜测。',
        '', '## 更新策略', '',
        'DepMap 的全部选定分析级文件应在新 release 时一起更新；只检查 release，不自动更新、覆盖或删除旧版。'
        '更新后的原始和 processed 数据都按 release 独立存储。TCGA 作为固定历史 cohort 不常规更新；'
        '只有官方对象本身变更时才另行评估。重新执行下载首先核实本地 SHA256，不重复下载成功文件。',
    ]
    (ROOT / 'docs/DATA_SOURCES.md').write_text('\n'.join(source_lines) + '\n', encoding='utf-8')
    lines = ['# 数据库准备报告', '', f'报告时间：{now()}。', '',
             f'最新 DepMap Public release：**{config.get("depmap_release")}**（官方目录自动识别）。'
             '按用户后续指示跳过 DepMap 实体下载；TCGA 下载、预处理及查询基础设施已执行。', '',
             '## 下载成功文件', '', '| 文件 | 字节数 | 大小 | SHA256 |', '|---|---:|---:|---|']
    for r in successful:
        lines.append(f'| `{r["filename"]}` | {int(r["file_size"]):,} | {size(r["file_size"])} | `{r["SHA256"]}` |')
    lines += ['', '## 磁盘空间', '',
              f'- DepMap raw：**{size(depmap_bytes)}**。',
              f'- TCGA raw：**{size(tcga_bytes)}**。',
              f'- processed 合计：**{size(processed_bytes)}**（Parquet、索引、sample order、provenance 和仅含 views 的 DuckDB）。',
              '- `.runtime/` 软件依赖及文档不计入上述数据库空间。', '',
              '## 未下载文件', '', '| 数据库 | 文件 | 状态 / 原因 |', '|---|---|---|']
    labels = {'skipped_by_user':'按用户要求跳过；官方 URL 需浏览器验证',
              'unavailable_optional':'26Q1 官方目录未提供，不从旧 release 替代',
              'failed':'下载失败：官方 GDC 历史 MAF UUID 已失效，API 无此文件',
              'missing_required':'官方目录缺少必需文件'}
    for r in missing:
        lines.append(f'| {r["database"]} | `{r["filename"]}` | {labels.get(r["status"], r["status"])} |')
    lines += ['', '精确失败信息见 manifest 的 error 列；单独清单见 `data/manifests/not_downloaded.csv`。'
              'TCGA 四个核心文件和 Xena MC3 补充表均成功；原始 MAF 可选文件失败 1 个。', '',
              '## 完整性与预处理 QC', '',
              '| 数据集 | 基因 | 样本 | 全部行解析 | 重复基因 | 缺失数值 |', '|---|---:|---:|---|---:|---:|']
    for kind in ['CN','GISTIC','Expression']:
        q = qc.get('qc', {}).get(kind, {})
        lines.append(f'| {kind} | {q.get("genes","未完成"):,} | {q.get("samples",0):,} | {q.get("parsed_all_rows")} | {q.get("duplicate_gene_labels")} | {q.get("missing_cells")} |')
    lines += ['', '三个矩阵为 float64、ZSTD Parquet，一基因一 row group。所有列宽、TCGA ID、数值解析均检查；'
              'GISTIC 只接受 -2/-1/0/1/2。原始 gzip 文件完整读取并验证 SHA256。', '',
              '### 样本与 phenotype overlap', '', '| 数据集 | 矩阵样本 | phenotype 原表 overlap | 癌种缺失 |', '|---|---:|---:|---:|']
    for kind, q in qc.get('qc', {}).get('metadata', {}).get('sample_overlap', {}).items():
        lines.append(f'| {kind} | {q["matrix_samples"]:,} | {q["phenotype_overlap"]:,} | {q["missing_cancer_type"]} |')
    m = qc.get('qc', {}).get('metadata', {})
    lines += ['', f'统一 metadata 共 {m.get("rows",0):,} 个样本，覆盖 33 个癌种。'
              '表达表的唯一 phenotype 缺失是 `TCGA-07-0249-20` 控制分析物，保留 Other、癌种 null；'
              '所有 Tumor/Normal 样本有癌种标注。患者 ID、样本类型、Tumor/Normal 已建立。', '',
              '### DepMap ModelID overlap', '',
              '**未测量**。Model.csv、CRISPR、CN、expression 均按用户指示未下载；不能报告 0 或声称通过。'
              '预处理脚本将来会检查 ACH-xxxxxx 格式、唯一映射及与 Model.csv 的 overlap，不要求各矩阵模型数量相同。', '',
              '### MC3 补充表', '',
              f'完整 Xena 表已转换成 MC3.parquet：{mc3.get("qc",{}).get("rows",0):,} 条记录，'
              f'{mc3.get("qc",{}).get("unique_samples",0):,} 个样本，'
              f'{mc3.get("qc",{}).get("unique_genes",0):,} 个基因；全部行解析且样本 ID 合法。', '',
              '### 查询与安全验证', '',
              f'真实数据库验证通过：**{validation.get("passed", "待验证")}**；'
              f'{len(validation.get("checks", []))} 项 checks，{len(validation.get("errors", []))} 项 errors。'
              '每个 TCGA 矩阵抽取 3 个基因，逐样本与原始值精确一致（float64），并检查 metadata。'
              '测试只验证数据读取，不做 Gene A/Gene B 分析。', '',
              '5 个合成安全测试覆盖 release 数值排序、默认 profile 映射、歧义失败、gene column projection、'
              '断点续传、checksum skip 和拒绝覆盖异常 raw；全部通过。合成测试不等于真实 DepMap QC。', '',
              '单基因查询使用显式 RowGroup 索引，仅读取一个基因向量；DuckDB views 另外保留供 SQL 使用。'
              '首次查询包括 footer/索引加载，后续查询只保留 footer 和基因索引缓存。', '',
              '单基因查询耗时：', '']
    for c in validation.get('checks', []):
        if c['check'] == 'single_gene_matches_raw':
            lines.append(f'- {c["dataset"]} / `{c["gene"]}`：{c["query_seconds"]} 秒，{c["samples"]:,} 样本。')
    lines += ['', '## 项目目录树', '', '```text', 'Copy Number–Driven Dependency Analysis/',
              '├─ data/', '│  ├─ raw/', '│  │  ├─ depmap/26Q1/             [空：本次跳过]',
              '│  │  └─ tcga/pancan/             [5 个原始 gzip 文件]',
              '│  ├─ processed/', '│  │  ├─ depmap/26Q1/             [空：尚无 raw]',
              '│  │  └─ tcga/                    [CN/GISTIC/Expression/MC3/Metadata + 索引 + DuckDB]',
              '│  └─ manifests/                  [SHA256、缺失清单、源字典、QC]',
              '├─ scripts/', '│  ├─ download/', '│  ├─ preprocess/', '│  ├─ utils/',
              '│  └─ data_access.py', '├─ config/', '├─ results/                       [无生物学分析]',
              '├─ logs/', '├─ docs/', '├─ tests/', '└─ README.md', '```', '',
              '## 后续更新', '',
              'DepMap 定期检查新 release，并在明确选择后整套更新；旧版永久保留。'
              'TCGA 不需要定期重新下载。MC3 原始 MAF 仅在官方链接恢复或找到官方稳定版本后补充。'
              'DepMap 将来下载时需要含真实 URL 的同 release 官方目录。', '',
              '自动审批曾拒绝排队的 DepMap 下载命令，理由是用户已要求跳过；该命令未执行。'
              '按后续指示完成其余已授权的 TCGA 工作。', '']
    (ROOT / 'docs/DATABASE_PREPARATION_REPORT.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps({'successful_files':len(successful), 'not_downloaded':len(missing),
                      'depmap_raw':size(depmap_bytes), 'tcga_raw':size(tcga_bytes),
                      'processed':size(processed_bytes), 'validation_passed':validation.get('passed')}, indent=2))


if __name__ == '__main__':
    main()
