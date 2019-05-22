## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  library(TCA)
#  library("data.table")
#  X <- data.frame(fread(file = "GSE42861_processed.subset.txt"), row.names=1)
#  W <- data.frame(fread(file = "GSE42861_processed.houseman_estimates.txt"), row.names=1)

## ----eval=FALSE----------------------------------------------------------
#  covars <- data.frame(
#    fread(file = "GSE42861_processed.covariates.txt"), row.names=1)
#  tech_var <- data.frame(
#    fread(file = "GSE42861_processed.control_probes_pcs.txt"),
#    row.names=1 )	# principal components calculated from control probes
#  C1 <- covars[,2:4]	# sex, age, and smoking status
#  C2 <- cbind(covars[,1],tech_var[,1:10])	# batch information and technical variation

## ----eval=FALSE----------------------------------------------------------
#  tca.mdl <- tca(X = X, W = W, C1 = C1, C2 = C2,
#                 parallel = TRUE, log_file = "GSE42861.tca.log")

## ----eval = FALSE--------------------------------------------------------
#  y <- as.matrix(covars[,5])

## ----eval = FALSE--------------------------------------------------------
#  C3 <- covars[,2:4]  # age, sex, and smoking status

## ----eval = FALSE--------------------------------------------------------
#  # if you wish to skip this step, simply load the pre-computed ReFACTor components
#  ref.scores <- data.frame(
#    fread(file = "GSE42861_processed.refactor_components.txt"), row.names=1)
#  
#  # otherwise, run the followings
#  
#  # load the full data matrix and run refactor
#  X.full <- data.frame(
#    fread(file = "GSE42861_processed.txt"), row.names=1)
#  ref <- refactor(X.full, k = 6, sparsity = 500,
#                  C = cbind(C1,C2), C.remove = TRUE,
#                  rand_svd = TRUE, log_file = "GSE42861.refactor.log")
#  ref.scores <- ref$scores

## ----eval = FALSE--------------------------------------------------------
#  C3 <- cbind(C3,ref$scores)

## ----eval = FALSE--------------------------------------------------------
#  cd4 <- colnames(W)[1] # the first column in W corresponds to CD4
#  res <- tcareg(X = X, tca.mdl = tca.mdl, y = y, C3 = C3,
#                test = "custom", null_model = NULL, alternative_model = c(cd4),
#                save_results = TRUE, output = "GSE42861.tcareg.CD4",
#                log_file = "GSE42861.tcareg.log",
#                features_metadata = "HumanMethylationSites.txt")

## ----eval = FALSE--------------------------------------------------------
#  plot(-log10(seq(from = 1/500, to = 1, by = 1/500)), -log10(sort(res$pvals)),
#       xlab = "Expected -log(p-val)", ylab = "Observed -log(p-val)")
#  abline(a=0, b=1); abline(h=-log10(0.05/500), col="red")

## ----eval = FALSE--------------------------------------------------------
#  hit.position <- order(res$pvals)[1] # 498
#  hit <- rownames(X)[hit.position] # cg11767757
#  hit.pval <- res$pvals[hit.position] # 7.95e-08

## ----eval = FALSE--------------------------------------------------------
#  # subset the estimates of the model parameters
#  mdl.tca.sub <- tcasub(tca.mdl, features = c(hit), log_file = "GSE42861.tcasub.log")
#  
#  # estimate cell-type-specific methylation for the associated site
#  Z_hat <- tensor(X = X[hit,], tca.mdl.sub, log_file = "GSE42861.tensor.log")

## ----eval = FALSE--------------------------------------------------------
#  cd4 <- Z_hat[[1]] # the first column in W is CD4, therefore take the first element in Z_hat
#  hit.cd4 <- cd4[1,]
#  d <- data.frame(hit.cd4, y)
#  boxplot(hit.cd4~y, data = d, main="CD4 methylation in cg11767757",
#          names=c("controls","cases"), ylab="Estimated methylation")

