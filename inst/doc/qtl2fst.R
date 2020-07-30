## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5)

## ---- load_packages-----------------------------------------------------------
library(qtl2)
library(qtl2fst)

## ---- load_iron_data----------------------------------------------------------
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

## ---- calc_alleleprob---------------------------------------------------------
pr <- calc_genoprob(iron, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)

## ----write_fst_db-------------------------------------------------------------
tmpdir <- file.path(tempdir(), "iron_genoprob")
dir.create(tmpdir)
fpr <- fst_genoprob(pr, "pr", tmpdir, quiet=TRUE)
fapr <- fst_genoprob(apr, "apr", tmpdir, quiet=TRUE)

## ----list_files---------------------------------------------------------------
list.files(tmpdir)

## ---- names_fpr---------------------------------------------------------------
names(fpr)

## ---- read_chromosome---------------------------------------------------------
apr_X <- fapr[["X"]]
dim(apr_X)

## -----------------------------------------------------------------------------
apr_X <- fapr$X
dim(apr_X)

## ----subset_fapr--------------------------------------------------------------
selected_ind <- subset(fapr, ind=1:20, chr=c("2","3"))
dim(fapr)

## ----subset_brackets----------------------------------------------------------
fapr_sub1 <- fapr[1:20, c("2","3")][["3"]]
fapr_sub2 <- fapr[,"2"]
fapr_sub23 <- fapr[,c("2","3")]
fapr_subX <- fapr[,"X"]

## ----select_markers-----------------------------------------------------------
dim(subset(fapr, mar=1:30))
dim(fapr[ , , dimnames(fapr)$mar$X[1:2]])

## ----cbind_fapr, warning=FALSE------------------------------------------------
fapr_sub223 <- cbind(fapr_sub2,fapr_sub23)

## ----rbind_fapr---------------------------------------------------------------
f23a <- fapr[1:20, c("2","3")]
f23b <- fapr[40:79, c("2","3")]
f23 <- rbind(f23a, f23b)

## ----subset_markers-----------------------------------------------------------
markers <- dimnames(fapr$X)[[3]][1:2]
dim(fapr[,,markers]$X)

## ----extract_markers_chr_X----------------------------------------------------
markers <- dimnames(fapr$X)[[3]]
dim(fapr$X[,,markers[1:2]])

## ----cbind_two_subsets--------------------------------------------------------
fapr2 <- fst_genoprob(subset(apr, chr="2"), "aprx", tmpdir, quiet=TRUE)
fapr3 <- fst_genoprob(subset(apr, chr="3"), "aprx", tmpdir, quiet=TRUE)
fapr32 <- cbind(fapr3,fapr2)
dim(fapr32)
list.files(tmpdir)

## ----names_fapr---------------------------------------------------------------
names(unclass(fapr))

## ----fapr_fst_path------------------------------------------------------------
unclass(fapr)$fst

## ----ind_chr_mar_pieces-------------------------------------------------------
sapply(unclass(fapr)[c("ind","chr","mar")], length)

## ----restore_from_subset------------------------------------------------------
fapr23 <- subset(fapr, chr=c("2","3"))
dim(fapr23)
dim(fst_restore(fapr23))

## ----path_to_fpr--------------------------------------------------------------
fst_path(fpr)

## ----new_path_to_fpr, warning=FALSE-------------------------------------------
fpr_newpath <- replace_path(fpr, tempdir())

## ----calc_genoprob_fst, warning=FALSE-----------------------------------------
fpr <- calc_genoprob_fst(iron, "pr", tmpdir, error_prob=0.002, overwrite=TRUE)

## ----genoprob_to_alleleprob_fst, warning=FALSE--------------------------------
fapr <- genoprob_to_alleleprob_fst(pr, "apr", tmpdir, overwrite=TRUE)

## ----genome_scan--------------------------------------------------------------
Xcovar <- get_x_covar(iron)
scan_pr <- scan1(fpr, iron$pheno, Xcovar=Xcovar)
find_peaks(scan_pr, iron$pmap, threshold=4)

## ----coef---------------------------------------------------------------------
coef16 <- scan1coef(fpr[,"16"], iron$pheno[,1])
blup16 <- scan1blup(fpr[,"16"], iron$pheno[,1])

## ----clean_up_files, include=FALSE--------------------------------------------
unlink(tmpdir)

