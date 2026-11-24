# env: env/pi/SCALPEL
library(data.table)
library(argparse)
library(Matrix)
parser <- ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('DGE')
parser$add_argument('PRED')
parser$add_argument('OUTP')
args <- parser$parse_args()
preds <- if (file.size(args$PRED) == 0) {
  data.table(bc=character(), gene_name=character(), transcript_name=character(), tr_prob=numeric())
} else {
  fread(args$PRED, header=FALSE, col.names=c('bc','gene_name','transcript_name','tr_prob'), colClasses=c('character','character','character','numeric'))
}
dge <- fread(args$DGE, check.names=FALSE)
features <- if (nrow(preds)) paste0(preds$gene_name, '***', preds$transcript_name) else character()
pairs <- unique(data.table(feature=features, gene=preds$gene_name, transcript=preds$transcript_name))
if (anyDuplicated(pairs$feature)) stop('Feature name collision: ', pairs$feature[anyDuplicated(pairs$feature)])
feature_names <- unique(features)
barcodes <- unique(preds$bc)
ri <- match(preds$gene_name, dge[[1L]])
ci <- match(preds$bc, names(dge)[-1L])
counts <- numeric(nrow(preds))
valid <- which(!is.na(ri) & !is.na(ci))
for (idx in split(valid, ci[valid])) {
  col <- ci[idx[1L]]
  counts[idx] <- dge[[col + 1L]][ri[idx]]
}
counts[is.na(counts)] <- 0
stem <- sub('\\.txt$', '', args$OUTP)
# An empty tr_prob is left by em_algorithm.R's left_join when it discards a gene ("Gene discarded", all fragment
# probabilities 0). Such isoforms get 0 and are recorded; NaN/Inf and other non-finite values are errors.
tr_prob <- preds$tr_prob
discarded <- is.na(tr_prob) & !is.nan(tr_prob)
tr_prob[discarded] <- 0
fwrite(data.table(bc=preds$bc[discarded], gene_name=preds$gene_name[discarded], transcript_name=preds$transcript_name[discarded], gene_count=counts[discarded]), paste0(stem, '_discarded.tsv'), sep='\t')
cat('EM-discarded records (empty tr_prob set to 0):', sum(discarded), '\n')
x <- tr_prob * counts
nonfinite <- which(!is.finite(x))
if (length(nonfinite)) stop('Non-finite tr_counts (', length(nonfinite), ' records), e.g. ', features[nonfinite[1L]], ' / ', preds$bc[nonfinite[1L]])
entries <- data.table(i=match(features, feature_names), j=match(preds$bc, barcodes), x=x)
conflicts <- entries[, .(values=uniqueN(x)), by=.(i,j)][values > 1L]
if (nrow(conflicts)) stop('Conflicting duplicate feature/barcode: ', feature_names[conflicts$i[1L]], ' / ', barcodes[conflicts$j[1L]])
entries <- unique(entries, by=c('i','j'))
mat <- sparseMatrix(i=entries$i, j=entries$j, x=entries$x, dims=c(length(feature_names), length(barcodes)), dimnames=list(feature_names, barcodes))
mat <- drop0(mat)
print(c(features=nrow(mat), barcodes=ncol(mat), nonzero=length(mat@x)))
writeMM(mat, args$OUTP)
writeLines(feature_names, paste0(stem, '_features.txt'))
writeLines(barcodes, paste0(stem, '_barcodes.txt'))
