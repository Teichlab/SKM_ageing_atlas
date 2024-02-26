if(1){
library(Matrix)
source("/nfs/team205/nk5/Team205/ed6/LMM.R")
# count 
Y=t(readMM("cpm.mtx.gz")) # original Y in mtx format
rownames(Y)=read.csv("features.csv",as.is=T)[,2]
#load("Y.Rbin") # R binary of Y above (transposed)

# meta data
mdata=read.csv("metadata.csv",as.is=T,header=T)
#names(mdata)[17]="celltype"
#mdata$percent_mito[mdata$percent_mito==0]=1e-5; mdata$percent_mito=log(mdata$percent_mito)
#mdata$percent_ribo[mdata$percent_ribo==0]=1e-5; mdata$percent_ribo=log(mdata$percent_ribo)
#mdata$n_counts=log(mdata$n_counts)
#mdata$n_genes=log(mdata$n_genes)

# cell proportion for each celltype
Z=model.matrix(~0+mdata$celltype)
P=(Y>0)%*%Z
P=t(t(P)/colSums(Z))
colnames(P)=gsub("mdata\\$celltype","",colnames(P))
save(P,file="P.Rbin")

# gene used (prop of cell expressing the gene > 5%)
flag=apply(P>0.05,1,sum)>0
Y = Y[flag,]
save(flag,file="flag_0.05_each_celltype.Rbin")
save(Y,P,Z,file="Y_0.05_each_celltype.Rbin")

# N genes and N UMIs
#mdata = cbind(mdata, n_gene=colSums(Y>0), n_umi=colSums(Y))

# linear mixed model
res=LFLMM(Y, mdata[,c("celltype", "X10X_version", "Ventilation_cat", "BMI", "DonorID","Age_bin", "Sex", "Donor_type")], ITRMAX=300, interactions=list(c("celltype","Age_bin")))
#Sample_ID is also a donorID

png(width=1000,height=1000,res=150,file="bar.png");par(mar=c(7,3,1,1));Barplot(res);dev.off()

save(res,file="res.Rbin")
}
# You need to parallelise below!!!!!
system("mkdir DE")
ucelltypes=unique(matrix(unlist(strsplit(colnames(res$Z)[res$H[,names(res$nh)=="celltype:Age_bin"]==1],":")),2)[1,])
for(i in ucelltypes){
	de = getBFInt(Y,res,"celltype:Age_bin",Celltype=i, DE1=NA)
	saveRDS(de, file=paste("DE/de_ctype2Age_group_",i,".RDS",sep=""))
}

# DE for marginal effects
de_celltype=getBF(Y,res,"celltype",DE1=1)
save(de_celltype,file="DE/de_celltype.Rbin")

## Collect and save DE results for each cell type
de_celltype_Age = (function(){
  gene=read.csv("features.csv",as.is=T)[flag,]
  mdata=read.csv("metadata.csv",as.is=T)
  celltypes=unique(mdata[,'celltype'])
  print(celltypes)
  celltype_list = flag5 = beta_old = beta_young = ltsr = gname = ensembl = nct = NULL
  th = 0.9
  for (celltype in celltypes){
    # DE result
    de=readRDS(paste("DE/de_ctype2Age_group_",celltype,".RDS",sep=""))
    # gene is expressed by at least 5% of cells
    if (ncol(de$beta)==2){
      flag5_1 = data.frame(flag5 = P[flag,celltype]>0.05)
      
      # DE genes i.e., Local True Sign Rate>th
      # flag_de = c(de$ltsr)>th & flag5
      # log Fold change
      #rownames(de$beta) = 1:nrow(de$beta) # replaced to be gene names
      beta_old1 = de$beta[,1,drop=F]
      beta_young1 = de$beta[,2,drop=F]
      # ltsr
      ltsr1 = de$ltsr[,1,drop=F]
      # collect
      colnames(beta_old1) = "beta_old"
      colnames(beta_young1) = "beta_young"
      colnames(ltsr1) = "ltsr"
      #gname = c(gname, gene[as.numeric(rownames(beta1)),2])
      gname = c(gname, gene[1:nrow(de$beta),"SYMBOL"])
      ensembl = c(ensembl, rownames(de$beta))
      beta_old = rbind(beta_old, beta_old1)
      beta_young = rbind(beta_young, beta_young1)
      ltsr = rbind(ltsr, ltsr1)
      flag5 = rbind(flag5, flag5_1)
      nct = c(nct, nrow(de$beta))
      celltype_list = c(celltype_list, celltype)
  }
}
data.frame(celltype=rep(celltype_list,nct),gene=gname, ENSEMBL = ensembl, ltsr=ltsr, beta_old=beta_old, beta_young=beta_young, flag5=flag5)
})()


write.table(de_celltype_Age,file="de_celltype_Age_interaction_all.txt",col=T,sep="\t",quote=F, row.names = FALSE)


