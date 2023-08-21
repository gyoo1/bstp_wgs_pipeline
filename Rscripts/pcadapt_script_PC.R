####################
# Run pcadapt scan #
####################
#
#       args[1]: zscores matrix
#       args[2]: prefix for output files
#

args = commandArgs(trailingOnly=TRUE)

library(RcppCNPy)
library(bigutilsr)


zscores <- Zscores
K <- ncol(zscores)

# For one component only
if (K == 1) {
        d2 <- (zscores - median(zscores))^2
} else {
        d2 <- dist_ogk(zscores)
}

write.table(pchisq(d2, df=K, lower.tail=F), file="PC.pcadapt.pval.txt", quote=F, row.names=F, col.names=F)
