# 本 fork 相对上游 SCALPEL 的改动

上游：<https://github.com/plasslab/SCALPEL>，基线提交 `3c31ffc`。本 fork（`yufan999/SCALPEL`）的 `main` 供 scAPA_Agent 的 `05_scalpel` 模块使用。以下只列相对 `3c31ffc` 的改动，可用 `git diff 3c31ffc main` 核对。

## 改动清单

| 文件 | 改动 |
|---|---|
| `workflows/reads_processing.nf` | Chromium 分支新增必填参数 `--dedup_mode`（`upstream` 或 `compound_cb_ub`）；samtools 的 index / sort / markdup（以及 `compound_cb_ub` 中的 view）使用 `task.cpus` 线程 |
| `src/combined_cell_umi_tag.py`（新增） | 给每条 read 生成临时标签 `XU:Z:<CB>_<UB>`，并统计缺少 CB 或 UB 的 reads，统计 JSON 输出到 stderr |
| `workflows/annotation_preprocessing.nf` | salmon 的 `-p` 由固定的 4 改为 `task.cpus` |
| `src/merge_DGE_PRED.R`、`src/APAtoseurat.R`、`workflows/isoform_quantification.nf` | DGE 改为稀疏矩阵，避免 `pivot_wider` 在大数据上整数溢出（见下） |
| `nextflow.config` | 资源设置改为作者机器的值：`params.cpus = 40`、`params.memory = 350.GB`、`time = '2000h'`、`small_rec` 的 `maxRetries = 1`（上游为 16 / 40 GB / 4h / 3） |

## dedup_mode

- `upstream`：与上游相同的去重命令 `samtools markdup -r --barcode-tag CB --barcode-tag UB`。上游重复写了 `--barcode-tag`，实际只有 UB 生效，这里保留不改。
- `compound_cb_ub`：先用 `src/combined_cell_umi_tag.py` 写入 `XU=CB_UB`，再 `markdup --barcode-tag XU`，最后用 `-x XU` 删除临时标签，使不同细胞中相同 UB 不会被互相当作重复。缺少 CB 或 UB 的 reads 会被丢弃并计数。
- 没有默认值：未传或取值非法会直接报错。
- 只改了 Chromium 分支，Drop-seq 分支仍是上游写法（同样重复设置 `--barcode-tag`）。
- 命令层面 `upstream` 与上游一致；两者的 BAM 级逐条对比尚未做。

## 保持与上游一致的部分

Chromium 分支始终使用 `-D CB:${bc_path}` 过滤细胞。未传 `--barcodes` 时，`bc_path` 是 cellranger 的 filtered barcodes。fork 历史里曾有一次提交（`ec87a64`）在未传 `--barcodes` 时去掉该过滤，现已改回上游行为。

## 稀疏 DGE

- 原因：`merge_DGE_PRED.R` 用 `pivot_wider` 生成稠密表，当细胞数乘特征数超过 R 的整数上限时会溢出。
- `<sample>_APADGE.txt` 改为 **MatrixMarket 稀疏矩阵**（文件名未改，但内容不再是制表符表格），并新增三个输出：
  - `<sample>_APADGE_features.txt`：行名（`gene***transcript`）
  - `<sample>_APADGE_barcodes.txt`：列名
  - `<sample>_APADGE_discarded.tsv`：被 EM 丢弃的记录（`bc`、`gene_name`、`transcript_name`、`gene_count`）
- EM 丢弃基因（`em_algorithm.R` 的 `left_join` 留下空 `tr_prob`）按 0 处理并写入 `_discarded.tsv`；`tr_prob * counts` 出现非有限值则报错退出。
- `APAtoseurat.R` 直接读稀疏矩阵构建 Seurat 对象。
- 下游若按制表符表格读取 `_APADGE.txt` 会出错，需要改读 MatrixMarket。

## 资源与运行时限

`nextflow.config` 里的 40 核、350 GB、2000 h 是作者机器的设置。其他机器应通过 `-c <file>` 覆盖；覆盖文件没有设置的项取 `nextflow.config` 的值。

## 在 scAPA_Agent 中的用法

scAPA_Agent 的 `modules/05_scalpel/05_1_SCALPEL_nextflow.sh` 以 `tools/SCALPEL/main.nf` 运行本仓库，用 `--dedup_mode` 传入去重模式，用 `-c config/05_scalpel/scalpel_resources.config` 覆盖资源。多文库联合运行（一张 samplesheet、单次 Nextflow 任务）在 Agent 一侧实现，不在本仓库。
