library(data.table)
library(Matrix)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE))
repo <- dirname(dirname(normalizePath(script)))
if (system2('git', c('-C', shQuote(repo), 'cat-file', '-e', shQuote('3c31ffc^{commit}')), stdout=FALSE, stderr=FALSE) != 0L) {
  cat('SKIP sparse DGE checks: upstream baseline 3c31ffc not available\n')
  quit(status=0)
}
rscript <- file.path(R.home('bin'), 'Rscript')
tmp <- tempfile('scalpel_sparse_')
dir.create(tmp)
original <- function(file) {
  dest <- file.path(tmp, basename(file))
  z <- system2('git', c('-C', shQuote(repo), 'show', shQuote(paste0('3c31ffc:', file))), stdout=dest)
  stopifnot(z == 0L)
  dest
}
oldmerge <- original('src/merge_DGE_PRED.R')
oldseurat <- original('src/APAtoseurat.R')
merge <- file.path(repo, 'src/merge_DGE_PRED.R')
seurat <- file.path(repo, 'src/APAtoseurat.R')
run <- function(script, args, success=TRUE, timed=FALSE) {
  timed <- timed && file.exists('/usr/bin/time')
  cmd <- if (timed) '/usr/bin/time' else rscript
  argv <- if (timed) c('-v', shQuote(rscript), shQuote(script), shQuote(args)) else c(shQuote(script), shQuote(args))
  status <- system2(cmd, argv)
  if (success && status != 0L || !success && status == 0L) stop('Unexpected exit status: ', status)
}
dge <- file.path(tmp, 'dge.txt'); pred <- file.path(tmp, 'pred.txt')
old <- file.path(tmp, 'old.txt'); out <- file.path(tmp, 'sparse.txt')
fwrite(data.table(GENE=c('gene_b','gene_a'), 'b-1'=c(3,4), 'a-1'=c(5,6)), dge, sep='\t')
p <- data.table(bc=c('b-1','a-1','missing-1','b-1','a-1','b-1'), gene_name=c('gene_b','gene_a','gene_b','missing','gene_b','gene_b'), transcript_name=c('t2','t1','t2','t3','t2','t2'), tr_prob=c(.5,.25,0,.3,0,.5))
writepred <- function(x) fwrite(x, pred, sep='\t', col.names=FALSE)
writepred(p)
run(oldmerge,c(dge,pred,old)); run(merge,c(dge,pred,out))
a <- as.matrix(read.table(old, header=TRUE, sep='\t', row.names=1, check.names=FALSE))
b <- as.matrix(readMM(out)); dimnames(b) <- list(readLines(file.path(tmp,'sparse_features.txt')),readLines(file.path(tmp,'sparse_barcodes.txt')))
stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=TRUE)))
# Both scripts write to the current directory, so run them in a temporary shell directory.
run_seurat <- function(script, input, sample) {
  cmd <- paste('cd', shQuote(tmp), '&&', shQuote(rscript), shQuote(script), shQuote(input), shQuote(sample), shQuote('unused'))
  stopifnot(system(cmd) == 0L)
  readRDS(file.path(tmp,paste0(sample,'_seurat.RDS')))
}
s1 <- run_seurat(oldseurat,old,'old'); s2 <- run_seurat(seurat,out,'sparse')
stopifnot(identical(colnames(s1),colnames(s2)), identical(rownames(s1),rownames(s2)), isTRUE(all.equal(as.matrix(s1[['RNA']]$counts),as.matrix(s2[['RNA']]$counts),tolerance=1e-12)))
cat('PASS small matrix and Seurat equivalence\n')
set.seed(1)
ng <- 60L; nc <- 80L; tpg <- 4L
genes <- sprintf('G_%02d', sample(ng)); cells <- sprintf('C%03d-1', sample(nc))
dm <- matrix(rpois(ng * nc, 2), ng, nc, dimnames=list(genes, cells))
dge2 <- file.path(tmp, 'dge2.txt'); pred2 <- file.path(tmp, 'pred2.txt')
fwrite(data.table(GENE=genes[-1L], dm[-1L, , drop=FALSE]), dge2, sep='\t')
grid <- CJ(bc=c(cells, 'extra-1'), g=genes, t=seq_len(tpg), sorted=FALSE)[sample(.N, 3000L)]
rp <- grid[, .(bc, gene_name=g, transcript_name=paste0('T', t), tr_prob=round(runif(.N), 3))]
fwrite(rp, pred2, sep='\t', col.names=FALSE)
run(oldmerge, c(dge2, pred2, old)); run(merge, c(dge2, pred2, out))
a <- as.matrix(read.table(old, header=TRUE, sep='\t', row.names=1, check.names=FALSE))
b <- as.matrix(readMM(out)); dimnames(b) <- list(readLines(file.path(tmp,'sparse_features.txt')), readLines(file.path(tmp,'sparse_barcodes.txt')))
stopifnot(isTRUE(all.equal(a, b, tolerance=1e-12, check.attributes=TRUE)), nrow(b) > 100L)
cat('PASS random medium matrix equivalence', dim(b), '\n')
# EM-discarded genes leave an empty tr_prob: the old flow kept NA (writeMM would turn NA into 1e308); expect 0 plus a record.
writepred(data.table(bc=c('b-1','b-1','a-1','a-1','a-1'), gene_name=c('gene_b','gene_a','gene_a','gene_a','gene_b'), transcript_name=c('t2','t1','t1','t5','t2'), tr_prob=c(.5,.25,NA,NA,.5)))
run(oldmerge,c(dge,pred,old)); run(merge,c(dge,pred,out))
a <- as.matrix(read.table(old, header=TRUE, sep='\t', row.names=1, check.names=FALSE))
b <- as.matrix(readMM(out)); dimnames(b) <- list(readLines(file.path(tmp,'sparse_features.txt')), readLines(file.path(tmp,'sparse_barcodes.txt')))
stopifnot(sum(is.na(a)) == 2L)
a[is.na(a)] <- 0
stopifnot(isTRUE(all.equal(a, b, tolerance=1e-12, check.attributes=TRUE)), all(is.finite(b)), max(abs(b)) < 1e300)
disc <- fread(file.path(tmp,'sparse_discarded.tsv'))
stopifnot(identical(disc$bc, c('a-1','a-1')), identical(disc$gene_name, c('gene_a','gene_a')), identical(disc$transcript_name, c('t1','t5')), all(disc$gene_count == 6))
s3 <- run_seurat(seurat,out,'discarded'); stopifnot(all(is.finite(s3[['RNA']]$counts@x)))
for (bad in c('NaN','Inf')) {
  writepred(data.table(bc=c('b-1','a-1'), gene_name=c('gene_b','gene_a'), transcript_name=c('t2','t1'), tr_prob=c(.5,0)))
  txt <- readLines(pred); txt[2] <- sub('\t0$', paste0('\t', bad), txt[2]); writeLines(txt, pred)
  run(merge,c(dge,pred,out),FALSE)
}
cat('PASS EM-discarded empty tr_prob -> 0 with record; NaN/Inf rejected\n')
writepred(rbind(p,p[1][, tr_prob := .75])); run(merge,c(dge,pred,out),FALSE)
writepred(data.table(bc=c('b-1','a-1'),gene_name=c('a***b','a'),transcript_name=c('c','b***c'),tr_prob=c(1,1))); run(merge,c(dge,pred,out),FALSE)
file.create(pred); writeLines(character(),pred); run(merge,c(dge,pred,out)); stopifnot(identical(dim(readMM(out)),c(0L,0L)), nrow(fread(file.path(tmp,'sparse_discarded.tsv')))==0L)
cat('PASS duplicates, collision, empty input\n')
big <- file.path(tmp, 'big')
dir.create(big, recursive=TRUE, showWarnings=FALSE)
tryCatch({
  nfeat <- 140000L; nbc <- 16600L
  bcs <- sprintf('BC%05d-1',seq_len(nbc))
  bigdge <- file.path(big,'dge.txt'); bigpred <- file.path(big,'pred.txt'); bigout <- file.path(big,'big_APADGE.txt')
  fwrite(as.data.table(c(list(GENE='G'),setNames(as.list(rep(1L,nbc)),bcs))),bigdge,sep='\t')
  fwrite(data.table(bc=rep(bcs,length.out=nfeat),gene_name='G',transcript_name=paste0('T',seq_len(nfeat)),tr_prob=.5),bigpred,sep='\t',col.names=FALSE)
  start <- proc.time()[['elapsed']]
  run(merge,c(bigdge,bigpred,bigout),timed=TRUE)
  timer <- if (file.exists('/usr/bin/time')) '/usr/bin/time -v' else ''
  cmd <- paste('cd',shQuote(big),'&&',timer,shQuote(rscript),shQuote(seurat),shQuote(bigout),'big unused')
  stopifnot(system(cmd)==0L)
  obj <- readRDS(file.path(big,'big_seurat.RDS'))
  stopifnot(identical(as.integer(dim(obj)),c(nfeat,nbc)), length(obj[['RNA']]$counts@x)==nfeat)
  cat('PASS big:',nfeat,'x',nbc,'nonzero',nfeat,'elapsed seconds',proc.time()[['elapsed']]-start,'\n')
}, finally=unlink(big,recursive=TRUE))
unlink(tmp,recursive=TRUE)
cat('PASS: all sparse DGE checks\n')
