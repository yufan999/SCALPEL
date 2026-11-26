# 本 fork 相对上游 SCALPEL 的改动

上游：<https://github.com/plasslab/SCALPEL>，基线提交 `3c31ffc`。本 fork（`yufan999/SCALPEL`）的 `main` 供 scAPA_Agent 的 `05_scalpel` 模块使用。以下只列相对 `3c31ffc` 的改动，可用 `git diff 3c31ffc main` 核对。

## 改动清单

| 文件 | 改动 |
|---|---|
| `workflows/reads_processing.nf` | Chromium 分支新增必填参数 `--dedup_mode`（`upstream` 或 `compound_cb_ub`）；samtools 的 index / sort / markdup（以及 `compound_cb_ub` 中的 view）使用 `task.cpus` 线程 |
| `src/combined_cell_umi_tag.py`（新增） | 给每条 read 生成临时标签 `XU:Z:<CB>_<UB>`，并统计缺少 CB 或 UB 的 reads，统计 JSON 输出到 stderr |
| `workflows/annotation_preprocessing.nf` | salmon 的 `-p` 由固定的 4 改为 `task.cpus` |
| `src/merge_DGE_PRED.R`、`src/APAtoseurat.R`、`workflows/isoform_quantification.nf` | DGE 改为稀疏矩阵，避免 `pivot_wider` 在大数据上整数溢出（见下）；`dge_generation` 的 `publishDir` 用 `saveAs` 把 `*_APADGE_discarded.tsv` 发布到 `<outputDir>/dge_discarded/`，其余输出仍在 `<outputDir>/` |
| `nextflow.config` | 资源设置改为作者机器的值：`params.cpus = 40`、`params.memory = 350.GB`、`time = '2000h'`、`small_rec` 的 `maxRetries = 1`（上游为 16 / 40 GB / 4h / 3） |

## 测试

- `Rscript tests/check_sparse_dge.R`：将稀疏 DGE 与上游 `3c31ffc` 的稠密脚本逐值对照，覆盖 EM-discarded、NaN/Inf、重复、空输入及 140000×16600 大矩阵。需要 Seurat、data.table、Matrix、argparse 等 R 包，且克隆中须有 `3c31ffc` 历史对象。
- `python3 tests/check_dedup_tagger.py`：用合成 BAM 检查 tagger 和 samtools markdup；PATH 中需要 samtools，bam2bed 可选。

## dedup_mode

- `upstream`：与上游相同的去重命令 `samtools markdup -r --barcode-tag CB --barcode-tag UB`。上游重复写了 `--barcode-tag`，实际只有 UB 生效，这里保留不改。
- `compound_cb_ub`：先用 `src/combined_cell_umi_tag.py` 写入 `XU=CB_UB`，再 `markdup --barcode-tag XU`，最后用 `-x XU` 删除临时标签，使不同细胞中相同 UB 不会被互相当作重复。缺少 CB 或 UB 的 reads 会被丢弃并计数。
- 没有默认值：未传或取值非法会直接报错。
- 只改了 Chromium 分支，Drop-seq 分支仍是上游写法（同样重复设置 `--barcode-tag`）。
- 2026-11-25 用 samtools 1.20 对比 HRR568828 的 chr21（464,015 条）和 chrM（6,073,212 条）：上游命令与 `upstream` 模式（`-@ 4`、`-@ 1` 各一次，再重跑一次）的 `samtools view` md5 和去掉 `@PG` 的 header md5 均相同，线程数 4 与 1 生成的 `.bai` 逐字节一致。`upstream` 与上游在 BAM 级逐条一致，线程参数不改变 `markdup` 保留的 reads。
- `compound_cb_ub` 与 `upstream` 有意不同：chr21 保留 463,777 条（`upstream` 为 464,015 条），chrM 保留 6,097,122 条（`upstream` 为 6,073,212 条）。

## 保持与上游一致的部分

Chromium 分支始终使用 `-D CB:${bc_path}` 过滤细胞。未传 `--barcodes` 时，`bc_path` 是 cellranger 的 filtered barcodes。fork 历史里曾有一次提交（`ec87a64`）在未传 `--barcodes` 时去掉该过滤，现已改回上游行为。

## 稀疏 DGE

- 原因：`merge_DGE_PRED.R` 用 `pivot_wider` 生成稠密表，当细胞数乘特征数超过 R 的整数上限时会溢出。
- `<sample>_APADGE.txt` 改为 **MatrixMarket 稀疏矩阵**（文件名未改，但内容不再是制表符表格），并新增三个输出：
  - `<sample>_APADGE_features.txt`：行名（`gene***transcript`）
  - `<sample>_APADGE_barcodes.txt`：列名
  - `<sample>_APADGE_discarded.tsv`：被 EM 丢弃的记录（`bc`、`gene_name`、`transcript_name`、`gene_count`），发布在 `<outputDir>/dge_discarded/` 下
- EM 丢弃的位置与上游不逐位等效：上游记为 `NA`（Seurat 计数矩阵中也是 `NA`），fork 记 0 并写入 `_discarded.tsv`。MatrixMarket 无法表示 `NA`，`Matrix::writeMM` 会将其写成 `1e308`；这是 fork 第一版稀疏化的缺陷，已在 `c230cad` 修复。`tr_prob * counts` 出现非有限值则报错退出。
- 与上游对比：2026-11-25 用 HRR568828、HRR568829 的同一份 `*.counts.txt` 与 `*_isoforms_quantified.txt` 分别运行上游稠密脚本和 fork 稀疏脚本，Seurat 计数矩阵维度（100018×1799、124814×2667）、行名列名及其顺序完全相同，稀疏结构相同，数值最大相对差为 2.2e-16。上游分别有 10 个、3 个 `NA`，位置与 `_discarded.tsv` 的 10 条、3 条逐一对应（原始 DGE 计数合计 5、4），fork 记 0。
- 影响：上游 `src/scalpel_library.R` 的 `Extract_aggTab` 先 `AggregateExpression` 再 `replace_na(., 0)`，一个 `NA` 细胞会让该文库该转录本的合计变 `NA` 再变 0，等于把整库计数清零；fork 记 0 不会。另外，下游若按上游方式保留 `NA`，Seurat 的 `CreateSeuratObject(min.features=…)` 会把含 `NA` 的细胞整个丢弃。
- `APAtoseurat.R` 直接读稀疏矩阵构建 Seurat 对象。
- 下游若按制表符表格读取 `_APADGE.txt` 会出错，需要改读 MatrixMarket。

## 资源与运行时限

`nextflow.config` 里的 40 核、350 GB、2000 h 是作者机器的设置。其他机器应通过 `-c <file>` 覆盖；覆盖文件没有设置的项取 `nextflow.config` 的值。

## 在 scAPA_Agent 中的用法

scAPA_Agent 的 `modules/05_scalpel/05_1_SCALPEL_nextflow.sh` 以 `tools/SCALPEL/main.nf` 运行本仓库，用 `--dedup_mode` 传入去重模式，用 `-c config/05_scalpel/scalpel_resources.config` 覆盖资源。多文库联合运行（一张 samplesheet、单次 Nextflow 任务）在 Agent 一侧实现，不在本仓库。

scAPA_Agent 的 `make test` 不再包含这两个测试，需要时在 `tools/SCALPEL` 下按上面的命令运行。

### 关于 DIU_table.csv

scAPA_Agent 运行不传 `--clusters`，`differential_isoform_usage` 因此令 `Clusters = orig.ident`，以每个文库为一组，合并文库内所有细胞后做卡方检验，输出 `DIU_table.csv` 和 `iDGE_seurat.RDS`。这是文库之间的探索性输出，没有样本层面的重复，p 值受细胞数影响，不用于病例对照结论；即使传入分组文件，仍是合并细胞后的检验。

scAPA_Agent 的正式样本层面分析在其 `05_3`：按（细胞类型 × 样本）做 pseudobulk，按供体配对，对供体做符号翻转置换检验；3 对样本只有 8 种符号翻转，p 值下限为 0.25，因此是描述性分析（候选基因由 TVD、dIR 效应量阈值定义）。
