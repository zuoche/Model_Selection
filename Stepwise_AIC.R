allele_subg <- function (dat, sSG, eSG, sSnp, eSnp, grnam, ...)
{
  snp<-c()
  SNP<-c()
  AICvalue <- matrix(ncol=1167, nrow=12)
  names(AICvalue) <- rep(c("KL_biogrp1","KL_biogrp2","KL_biogrp3","KL_biogrp4","KL_biogrp5","KL_biogrp6","KL_biogrp7","KL_biogrp8","KL_biogrp9","KL_biogrp10","KL_biogrp11","KL_biogrp12"), 1167)
  
  da<-dat[,c(sSG:eSG,sSnp:eSnp)]
  lend<-eSG-sSG+1
  lend2<-eSnp-sSnp+1
  
  for (i in (lend+1):(lend2+lend)) {
    ssnp<-names(da)[i]
    formul<-reformulate(termlabels = grnam, response = ssnp)
    mod <-multinom(formul, data=da, model=T)
    snp<-append(snp,ssnp,after=length(snp))
    
    aic <- stepAIC(mod, direction="backward")
    
    for (j in (i*12-155):(i*12-144)) ifelse(names(AICvalue)[j] %in% summary(aic)$coefnames, AICvalue[j] <- summary(aic)$AIC, AICvalue[j] <- NA)
    
  }
  
  for (i in 1:length(snp)) SNP[i] <- substr(snp[i], 1, nchar(as.character(snp[i]))-2)
  
  
  return(AICvalue)
  
}

test<-allele_subg(ccp2,3,14,21,1187,"KL_biogrp1+KL_biogrp2+KL_biogrp3+KL_biogrp4+KL_biogrp5+KL_biogrp6+KL_biogrp7+KL_biogrp8+KL_biogrp9+KL_biogrp10+KL_biogrp11+KL_biogrp12")
