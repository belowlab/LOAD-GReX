#! /usr/bin/evn Rscript
args= commandArgs(trailingOnly=T)
cohort=args[1]
print(paste("working on ",cohort,"...",sep=""))
library("logistf")
library("data.table")
library("R.utils")
tissue=read.table(args[2], header=F)
cov = read.table(paste('/data30t1/LOAD/predixcan_ADGC/cov/',cohort,'/predixcan_cov.txt',sep=""),header=T,stringsAsFactors=F)
cov$Pheno= as.integer(as.character(cov$Pheno))-1
fam = read.table(paste('/data30t1/LOAD/predixcan_ADGC/cov/',cohort,'/fam.txt',sep=""),header=F, col.names=c("FID","IID","PID","MID","SEX","PHENO"), stringsAsFactors=F)
tiss = read.table('mashr_db.txt',header=F)
cond_ens = read.table('/data30t1/LOAD/predixcan_ADGC/cond/cond_gene_grch38.txt', header=T, sep = '\t', stringsAsFactors=F)
gene_list = read.table('/data30t1/LOAD/predixcan_ADGC/cond/gene_for_adj_grch38.txt', header=T, sep='\t', stringsAsFactors = F)
cond_snp = read.table(paste('/data30t1/LOAD/predixcan_ADGC/tagSNP/data_snp/',cohort,'.dosage', sep=""), header=F)
snp = cbind(fam[,c('FID','IID')], t(cond_snp[,c(7:length(cond_snp[1,]))]))
colnames(snp) = c('FID','IID',as.character(cond_snp$V2))
tagsnp = read.table(paste('/data30t1/LOAD/predixcan_ADGC/tagSNP/',cohort,'/TagSNP.txt', sep=""), header=F, col.names=c('target','tag','r2','pos'))
tagsnp$tag = as.character(tagsnp$tag)
tagsnp$target = as.character(tagsnp$target)
firth_cond_snp = function(ens, chr, start, end){
        print(ens)
        start = as.numeric(start)
        end = as.numeric(end)
        data=merge(predix[,c('FID','IID',genes[ens])], cov, by=c('FID','IID'), all.y=T)
        if(sd(data[which(data$Pheno==0),genes[ens]])==0 || sd(data[which(data$Pheno==1),genes[ens]])==0 || chr==16){
        c(Ensembl.ID=genes[ens], effect_size=NA, pvalue=NA, n=NA, Ref=NA, nonRef=NA, adj=NA)
        }else {
        for(tsnp in unique(as.character(cond_ens$SNP[which(cond_ens$chr==chr)]))){
                if((start-1e7) < cond_ens$pos[cond_ens$SNP==tsnp][1] && (end+1e7)> cond_ens$pos[cond_ens$SNP==tsnp][1]){
                        if(tagsnp$tag[which(tagsnp$target==tsnp)] %in% colnames(snp)){
                        data = merge(data, snp[,c('FID','IID', tagsnp$tag[which(tagsnp$target==tsnp)])], by=c("FID","IID"))}}}
        form = as.formula(paste('Pheno ~ ', paste(setdiff(colnames(data),c('FID','IID','Pheno','pc4','omit')), collapse =' + ', sep=' '), sep=''))
        firth =logistf(formula = form, family = "binomial", data = data)
        adjtest = setdiff(colnames(data),c('FID','IID','Pheno','Age','Sex','pc1','pc2','pc3','pc4','omit',genes[ens]))
        adj = paste(setdiff(colnames(data),c('FID','IID','Pheno','Age','Sex','pc1','pc2','pc3','pc4','omit',genes[ens])), collapse =',')
        if(length(adjtest)==0){
        c(Ensembl.ID=genes[ens], effect_size=NA, pvalue=NA, n=NA, Ref=NA, nonRef=NA, adj=NA)
        } else {
        c(Ensembl.ID=genes[ens] ,effect_size=firth$coefficients[2], pvalue=firth$prob[2]+1e-15, n= length(firth$weights),Ref='A', nonRef="T", adj= adj)
}}}
firth_cond_geno = function(ens, chr, start, end){
        print(ens)
        start = as.numeric(start)
        end = as.numeric(end)
        data=merge(predix[,c('FID','IID',genes[ens])], cov, by=c('FID','IID'), all.y=T)
        if(sd(data[which(data$Pheno==0),genes[ens]])==0 || sd(data[which(data$Pheno==1),genes[ens]])==0){
        c(Ensembl.ID=genes[ens], effect_size=NA, pvalue=NA, n=NA, Ref=NA, nonRef=NA, adj= NA)
        }else {
        for(adgene in unique(as.character(cond_ens$ens[which(cond_ens$chr==chr)]))){
                if((start-1e7) < max(cond_ens$pos[cond_ens$ens==adgene]) && (end+1e7)> min(cond_ens$pos[cond_ens$ens==adgene])){
                        if(genes[adgene] %in% colnames(alladgene) && sd(alladgene[,genes[adgene]])>0 ){
                        data = merge(data, alladgene[,c('FID','IID', genes[adgene])], by=c("FID","IID"))}}}
        form = as.formula(paste('Pheno ~ ', paste(setdiff(colnames(data),c('FID','IID','Pheno','pc4','omit')), collapse =' + ', sep=' '), sep=''))
        firth =logistf(formula = form, family = "binomial", data = data)
        adjtest = setdiff(colnames(data),c('FID','IID','Pheno','Age','Sex','pc1','pc2','pc3','pc4','omit',genes[ens]))
        adj = paste(setdiff(colnames(data),c('FID','IID','Pheno','Age','Sex','pc1','pc2','pc3','pc4','omit',genes[ens])), collapse =',')
        if(length(adjtest)==0){
        c(Ensembl.ID=genes[ens], effect_size=NA, pvalue=NA, n=NA, Ref=NA, nonRef=NA, adj= NA)
        } else {
        c(Ensembl.ID=genes[ens] ,effect_size=firth$coefficients[2], pvalue=firth$prob[2]+1e-22, n= length(firth$weights), Ref='A', nonRef="T", adj= adj)
}}}
firth = function(ens){
        print(ens)
        data=merge(predix[,c('FID','IID', ens)], cov, by=c('FID','IID'), all.y=T)
        if(sd(data[which(data$Pheno==0), ens])==0 || sd(data[which(data$Pheno==1), ens])==0){
        c(Ensembl.ID=ens, effect_size=NA, pvalue=NA, n=NA, Ref=NA, nonRef=NA)
        }else {
        form = as.formula(paste('Pheno ~ ', paste(setdiff(colnames(data),c('FID','IID','Pheno','pc4','omit')), collapse =' + ', sep=' '), sep=''))
        firth =logistf(formula = form, family = "binomial", data = data)
        c(Ensembl.ID=ens ,effect_size=firth$coefficients[2], pvalue=firth$prob[2]+1e-22, n= length(firth$weights), Ref='A', nonRef="T")}
}
for(i in c(1:length(tissue[,1]))){
        print(paste("working on ",cohort," with ",tissue[i,], sep=""))
        predix = fread(paste("/data30t1/LOAD/predixcan_ADGC/gtex.v8/predixcan/",cohort,"/",tissue[i,],".txt.gz",sep=""),header=T, data.table=F)
        genes = colnames(predix)[c(-1, -2)]
#       predix = cbind(fam[,c('FID','IID')], predix)
        result = t(sapply(genes, firth))
        colnames(result)=c("Ensembl.ID","effect_size","pvalue","N","Ref","nonRef")
        tt = tissue[i,]
#       tt=gsub("_imputed_europeans_tw_0.5_signif.db","",tissue[i,])
#       tt=gsub("_0.5.db","",tt)
#       tt=gsub("_newMetax.db","",tt)
        write.table(result, paste("./uncond/",tt,"/",cohort,".firth.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
        print('finished unconditional analysis')
        names(genes)=sapply(genes, function(x) strsplit(x, split='.',fixed = TRUE)[[1]][1])
        candi.gene = intersect(names(genes), gene_list$gene)
        predix =  predix[,c('FID','IID',genes[candi.gene])]
        result = t(sapply(candi.gene, function(x) firth_cond_snp(x, gene_list[which(gene_list==x),'chr'],gene_list[which(gene_list==x),'start'], gene_list[which(gene_list==x),'end'])))
        colnames(result)=c("Ensembl.ID","effect_size","pvalue","N","Ref","nonch(gene_list==x),'end'])))Ref","adj_SNP")
        write.table(result, paste("./cond_snp/",tt,"/",cohort,"_firth.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
        print('finished conditional for SNP')
        alladgene= cbind(fam[,c('FID','IID')], predix[,genes[intersect(cond_ens$ens,names(genes))]])
#       print(colnames(alladgene))
        predix = predix[,c('FID','IID',setdiff(colnames(predix), colnames(alladgene)))]
#       genes = colnames(predix)[c(-1,-2)]
#       names(genes)=sapply(genes, function(x) strsplit(x, split='.',fixed = TRUE)[[1]][1])
        candi.gene = colnames(predix)[c(-1,-2)]
        candi.gene = sapply(candi.gene, function(x) strsplit(x, split='.',fixed = TRUE)[[1]][1])
        result = t(sapply(candi.gene, function(x) firth_cond_geno(x, gene_list[which(gene_list==x),'chr'],gene_list[which(gene_list==x),'start'], gene_list[which(gene_list==x),'end'])))
        colnames(result)=c("Ensembl.ID","effect_size","pvalue","N","Ref","nonRef","adj_gene")
        write.table(result, paste("./cond_gene/",tt,"/",cohort,"_firth.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
        print('finished conditional for gene')
}
