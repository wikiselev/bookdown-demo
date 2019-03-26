run_MetaNeighbor_US<-function(vargenes, data, celltypes, pheno){
  
  cell.labels=matrix(0,ncol=length(celltypes),nrow=dim(pheno)[1])
  rownames(cell.labels)=colnames(data)
  colnames(cell.labels)=celltypes
  for(i in 1:length(celltypes)){
    type=celltypes[i]
    m<-match(pheno$Celltype,type)
    cell.labels[!is.na(m),i]=1
  }
  
  m<-match(rownames(data),vargenes)
  cor.dat=cor(data[!is.na(m),],method="s")
  rank.dat=cor.dat*0
  rank.dat[]=rank(cor.dat,ties.method="average",na.last = "keep")
  rank.dat[is.na(rank.dat)]=0
  rank.dat=rank.dat/max(rank.dat)
  sumin    =  (rank.dat) %*% cell.labels
  sumall   = matrix(apply(rank.dat,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])
  predicts = sumin/sumall
  
  cell.NV=matrix(0,ncol=length(celltypes),nrow=length(celltypes))
  colnames(cell.NV)=colnames(cell.labels)
  rownames(cell.NV)=colnames(cell.labels)
  
  for(i in 1:dim(cell.labels)[2]){
    predicts.temp=predicts
    m<-match(pheno$Celltype,colnames(cell.labels)[i])
    study=unique(pheno[!is.na(m),"Study_ID"])
    m<-match(pheno$Study_ID,study)
    pheno2=pheno[!is.na(m),]
    predicts.temp=predicts.temp[!is.na(m),]
    predicts.temp=apply(abs(predicts.temp), 2, rank,na.last="keep",ties.method="average")
    filter=matrix(0,ncol=length(celltypes),nrow=dim(pheno2)[1])
    m<-match(pheno2$Celltype,colnames(cell.labels)[i])
    filter[!is.na(m),1:length(celltypes)]=1
    negatives = which(filter == 0, arr.ind=T)
    positives = which(filter == 1, arr.ind=T)
    predicts.temp[negatives] <- 0
    np = colSums(filter,na.rm=T)
    nn = apply(filter,2,function(x) sum(x==0,na.rm=T))
    p =  apply(predicts.temp,2,sum,na.rm=T)
    cell.NV[i,]= (p/np - (np+1)/2)/nn
  }
  
  cell.NV=(cell.NV+t(cell.NV))/2
  return(cell.NV)
  
}

get_variable_genes<-function(data, pheno) {
var.genes1=vector("list")
experiment=unique(pheno$Study_ID)
j=1
for(exp in experiment){
  dat.sub=data[,pheno$Study_ID==exp]
    genes.list=vector("list")
    med.dat=apply(dat.sub,1,median)
    var.dat=apply(dat.sub,1,var)
    quant.med=unique(quantile(med.dat,prob=seq(0,1,length=11),type=5))
    genes.list=vector("list",length=length(quant.med))
    for(i in 1:length(quant.med)){
      if(i==1){
        filt1=med.dat<=quant.med[i]
        var.temp=var.dat[filt1]
        quant.var=quantile(var.temp,na.rm=T)
        filt2=var.temp>quant.var[4]
        genes.list[[i]]=names(var.temp)[filt2]
      }
      else {
        filt1=med.dat<=quant.med[i]&med.dat>quant.med[i-1]
        var.temp=var.dat[filt1]
        quant.var=quantile(var.temp,na.rm=T)
        filt2=var.temp>quant.var[4]
        genes.list[[i]]=names(var.temp)[filt2]
      }
    }
    temp=length(genes.list)
    var.genes1[[j]]=unlist(genes.list[1:temp-1])
    j=j+1
}
var.genes=Reduce(intersect, var.genes1)
return(var.genes)
}


get_top_hits <- function(cell.NV, pheno, threshold=0.95, filename) {
  
  type_by_study=table(pheno[,c("Celltype","Study_ID")])
  m<-match(rownames(cell.NV),rownames(type_by_study))
  f.a=!is.na(m)
  f.b=m[f.a]
  cell.NV=cell.NV[f.a,f.a]
  type_by_study=type_by_study[f.b,]
  
  for(i in 1:dim(type_by_study)[2]){
    filt=type_by_study[,i]!=0
    cell.NV[filt,filt]=0
  }
  
  diag(cell.NV)=0
  temp=vector()
  for(i in 1:dim(cell.NV)[1]){
    temp=c(temp,which.max(cell.NV[i,]))
  }
  temp=cbind(rownames(cell.NV),temp)
  for(i in 1:dim(cell.NV)[1]){
    temp[i,2]=cell.NV[i,as.numeric(temp[i,2])]
  }
  
  recip=temp[duplicated(temp[,2]),]
  filt=as.numeric(temp[,2])>=threshold
  recip=rbind(recip,temp[filt,])
  recip=cbind(recip,c(rep("Reciprocal_top_hit",each=dim(recip)[1]-sum(filt)),rep(paste("Above",threshold,sep="_"),each=sum(filt))))
  recip=recip[!duplicated(recip[,2]),]
  
  recip2=cbind(rownames(recip),recip[,1:3])
  colnames(recip2)=c("Celltype_1","Celltype_2","Mean_AUROC","Match_type")
  rownames(recip2)=NULL
  recip=recip2[order(recip2[,3],decreasing=T),]
  recip2=as.data.frame(recip)
  recip2[,3]=round(as.numeric(as.character(recip2[,3])),2)
  write.table(recip,file=filename,sep="\t",quote=F)
  return(recip2)
}


