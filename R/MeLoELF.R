### GENERAL INFO
#
# An R script for analyzing MeLoELF-seq data (Nanopore sequencing of oligomerized methylation reaction products)
#
# Parsing/alignment by Wayne Hemphill
# Processivity analysis by Ella R. Tommer
#
### COMMENTS
#
# On a 16-core Apple M3 Max, maps about 9,000 reads/min with full parallelization
#
###############

### Script Parameterization (edit as necessary)
MeLoELF <- function(parent,
                    target,
                    FWD.sites,
                    REV.sites,
                    mdir=getwd(),
                    sam.file='auto',
                    sam.indices=c(10,33,34,0),
                    crunch.too=T,
                    process=T,
                    pre.ligated=F,
                    methyl.type='B',
                    thresh.meth='BM',
                    BM.strand.data='r',
                    var.thresh=0,
                    T1.err=0.01,
                    T2.err=0.05,
                    target_fwd_auc=0.9,
                    read.length=NULL,
                    padding=2,
                    completeness=0.7,
                    matching=0.9,
                    ins.tol=1,
                    modk.fillC=F,
                    exact.search=F,
                    align.file='align.RData',
                    processed.file='processed.RData',
                    seq.file='sequences.txt',
                    melo.file='methlocations.txt',
                    meca.file='methcalls.txt',
                    plot_title='',
                    plot.nom=c('FWD','REV')){


#######################################################################
### Script (don't edit past here)
#######################################################################

setwd(mdir)

#######################
## Give errors & warnings for input
#######################

# ERRORS
if(nchar(parent)!=nchar(target)){
  stop('ERROR! The inputted FWD/parent and REV/target reference sequences are different lengths. MeLoELF only supports analysis of blunt-end substrates.')
}

# WARNINGS
if(exact.search==T){
  warning('Motif-based search enabled...the completeness, matching, and padding parameters will be ignored.')
}
if(pre.ligated){
  warning('Fragment re-oligomerization enabled...the completeness and matching parameters will be repurposed.')
}
if(!is.null(read.length) & !(exact.search==T)){
  warning('WARNING:  Mapping long reads without the motif search algorithm drastically increases computational time.')
}
if(var.thresh==2){
  warning('Per-site thresholding enabled...the BM.strand.data parameter will be ignored.')
}
if(length(FWD.sites) != length(REV.sites)){
  warning('WARNING:  The FWD/parent and REV/target strands have different numbers of methylation sites...some analyses may be affected adversely.')
}
if(Sys.info()['sysname']=='Windows'){
  warning('WARNING! -- The MeLoELF package was designed with Unix systems in mind, and may not work on Windows PCs. The sam -> txt file processing invokes bash/awk functionality, so at minimum, the sequence, methlocations, and methcalls txt files must be supplied manually to the mdir directory.')
}
if(grepl('fastq|fq',sam.file,ignore.case = T)){
  warning('FASTQ file provided...it will be treated as bisulfite sequencing data, and input adapted accordingly. Adapting (~50k reads/min)...')
  used.fastq=T
  if(!modk.fillC){
    stop("ERROR: Bisulfite sequencing FASTQ files should be analyzed with modk.fillC=T to ensure proper methylation assignment for C basecalls.")
  }
}else{
  used.fastq=F
}
  
#######################
## Load necessary non-base packages
#######################

library(stringr)
library(foreach)
library(doParallel)
library(pracma)
library(data.table)
library(dplyr)
library(MASS)

#######################
## Create custom functions for later analysis
#######################

# Indel repair algorithm
indel.fix <- function(DATA,FWD,REV){
  find.indel <- function(FWD.I,REV.I,fID){
    res=data.frame('indel'=rep(NA,times=length(fID)),'id'=rep(NA,times=length(fID)))
    for(i in 1:(length(fID)-1)){
      if(!is.na(fID[i]) & !is.na(fID[i+1])){
        if(fID[i]==fID[i+1]){
          if(sum(which(!is.na(c(FWD.I[i,],REV.I[i,]))) %in% which(!is.na(c(FWD.I[i+1,],REV.I[i+1,]))))==0){
            if(mean(which(!is.na(c(FWD.I[i,],REV.I[i,])))) < mean(which(!is.na(c(FWD.I[i+1,],REV.I[i+1,]))))){
              if(diff(range(na.omit(c(FWD.I[i:(i+1),],REV.I[i:(i+1),]))))<=(ncol(FWD.I)+ins.tol)){
                res$indel[i]=T
                res$id[i]=fID[i]
              }else{
                res$indel[i]=F
              }
            }else{
              res$indel[i]=F
            }
          }else{
            res$indel[i]=F
          }
        }else{
          res$indel[i]=F
        }
      }else{
        res$indel[i]=F
      }
    }
    return(res)
  }

  sort.id=order(rowMeans(cbind(DATA$FWD.I,DATA$REV.I),na.rm = T))
  #
  FWD.index=DATA$FWD.I[sort.id,]
  REV.index=DATA$REV.I[sort.id,]
  Q.reads=DATA$Q[sort.id,]
  #
  indels=find.indel(FWD.index,REV.index,Q.reads[,3])
  indel.id=which(indels$indel)
  if(isempty(indel.id)){
    return(DATA)
  }else{
    FWD.align=DATA$FWD.align[sort.id,]
    REV.align=DATA$REV.align[sort.id,]
    FWD.Chm=DATA$FWD.Chm[sort.id,]
    REV.Chm=DATA$REV.Chm[sort.id,]
    FWD.Cm=DATA$FWD.Cm[sort.id,]
    REV.Cm=DATA$REV.Cm[sort.id,]
    #
    for(i in length(indel.id):1){
      if((indels$id[which(indels$indel)])[i]=='FWD'){
        FWD.align[indel.id[i],which(FWD.align[indel.id[i]+1,]!='x')]=FWD.align[indel.id[i]+1,which(FWD.align[indel.id[i]+1,]!='x')]
        FWD.index[indel.id[i],which(!is.na(FWD.index[indel.id[i]+1,]))]=FWD.index[indel.id[i]+1,which(!is.na(FWD.index[indel.id[i]+1,]))]
        FWD.Chm[indel.id[i],which(!is.na(FWD.Chm[indel.id[i]+1,]))]=FWD.Chm[indel.id[i]+1,which(!is.na(FWD.Chm[indel.id[i]+1,]))]
        FWD.Cm[indel.id[i],which(!is.na(FWD.Cm[indel.id[i]+1,]))]=FWD.Cm[indel.id[i]+1,which(!is.na(FWD.Cm[indel.id[i]+1,]))]
        Q.reads[indel.id[i],1]=(diff(range(which(FWD.align[indel.id[i],]!='x')))+1)/length(FWD)
        Q.reads[indel.id[i],2]=sum(FWD.align[indel.id[i],]==FWD)/(diff(range(which(FWD.align[indel.id[i],]!='x')))+1)
      }
      if((indels$id[which(indels$indel)])[i]=='REV'){
        REV.align[indel.id[i],which(REV.align[indel.id[i]+1,]!='x')]=REV.align[indel.id[i]+1,which(REV.align[indel.id[i]+1,]!='x')]
        REV.index[indel.id[i],which(!is.na(REV.index[indel.id[i]+1,]))]=REV.index[indel.id[i]+1,which(!is.na(REV.index[indel.id[i]+1,]))]
        REV.Chm[indel.id[i],which(!is.na(REV.Chm[indel.id[i]+1,]))]=REV.Chm[indel.id[i]+1,which(!is.na(REV.Chm[indel.id[i]+1,]))]
        REV.Cm[indel.id[i],which(!is.na(REV.Cm[indel.id[i]+1,]))]=REV.Cm[indel.id[i]+1,which(!is.na(REV.Cm[indel.id[i]+1,]))]
        Q.reads[indel.id[i],1]=(diff(range(which(REV.align[indel.id[i],]!='x')))+1)/length(REV)
        Q.reads[indel.id[i],2]=sum(REV.align[indel.id[i],]==REV)/(diff(range(which(REV.align[indel.id[i],]!='x')))+1)
      }
    }
    #
    FWD.align[(indel.id+1),]='x'
    REV.align[(indel.id+1),]='x'
    FWD.Chm[(indel.id+1),]=NA
    REV.Chm[(indel.id+1),]=NA
    FWD.Cm[(indel.id+1),]=NA
    REV.Cm[(indel.id+1),]=NA
    FWD.index[(indel.id+1),]=NA
    REV.index[(indel.id+1),]=NA
    Q.reads[(indel.id+1),]=NA
    #
    res=list('FWD.align'=FWD.align,'REV.align'=REV.align,'FWD.Chm'=FWD.Chm,'FWD.Cm'=FWD.Cm,'REV.Chm'=REV.Chm,'REV.Cm'=REV.Cm,'Q'=as.matrix(Q.reads),'FWD.I'=FWD.index,'REV.I'=REV.index)
    return(res)
  }
}

# Exact motif search algorithm
find.motif <- function(read,Cm,Chm,C.key,read.length,FWD,REV){

  if(nchar(read)<min(read.length) | nchar(read)>max(read.length)){
    DATA='blank' # bypasses polymers outside desired length range
  } else {
    if((length(Cm)+length(Chm))!=length(C.key)){
      warning('WARNING! Read MM and ML tag info are different lengths...skipping read!')
      DATA='blank'
    }else{
      polymer=str_split(read,pattern = '')[[1]] # takes the polymer sequence and converts it from a single string into a vector of 1 base per index
      Cm.ids=which(polymer=='C')[Cm]
      Cm.scores=C.key[(length(Chm)+1):(length(Chm)+length(Cm))]
      Chm.ids=which(polymer=='C')[Chm]
      Chm.scores=C.key[1:length(Chm)]
      #
      FWD.match=which(polymer==FWD[1])
      REV.match=which(polymer==REV[1])
      for(i in 2:length(FWD)){
        FWD.match=FWD.match[which(polymer[FWD.match+i-1] == FWD[i])]
        REV.match=REV.match[which(polymer[REV.match+i-1] == REV[i])]
      }
      #
      if(length(c(FWD.match,REV.match))>0){
        FWD.align=matrix(c(rep(FWD,times=length(FWD.match)),rep('x',times=length(FWD)*length(REV.match))),ncol = length(FWD),nrow = length(c(FWD.match,REV.match)),byrow = T)
        colnames(FWD.align) <- paste0(FWD,'.',1:length(FWD))
        REV.align=matrix(c(rep('x',times=length(REV)*length(FWD.match)),rep(REV,times=length(REV.match))),ncol = length(REV),nrow = length(c(FWD.match,REV.match)),byrow = T)
        colnames(REV.align) <- paste0(REV,'.',1:length(REV))
        #
        quality=cbind('comp'=rep(1,times=length(c(FWD.match,REV.match))),'match'=rep(1,times=length(c(FWD.match,REV.match))),'id'=rep(c('FWD','REV'),times=c(length(FWD.match),length(REV.match))))
        #
        FWD.I=matrix(c(FWD.match-1,rep(NA,times=length(REV.match))),nrow = length(c(FWD.match,REV.match)),ncol = length(FWD))+matrix(c(rep(1:length(FWD),times=length(FWD.match)),rep(0,times=length(REV.match)*length(FWD))),nrow = length(c(FWD.match,REV.match)),ncol = length(FWD),byrow = T)
        REV.I=matrix(c(rep(NA,times=length(FWD.match)),REV.match-1),nrow = length(c(FWD.match,REV.match)),ncol = length(REV))+matrix(c(rep(0,times=length(FWD.match)*length(REV)),rep(1:length(REV),times=length(REV.match))),nrow = length(c(FWD.match,REV.match)),ncol = length(REV),byrow = T)
        #
        FWD.Chm=t(matrix(sqrt((FWD=='C')-1),ncol = length(FWD),nrow = length(c(FWD.match,REV.match)),byrow = T))
        REV.Chm=t(matrix(sqrt((REV=='C')-1),ncol = length(REV),nrow = length(c(FWD.match,REV.match)),byrow = T))
        if(!modk.fillC){
          FWD.Chm[FWD.Chm==0]=NA
          REV.Chm[REV.Chm==0]=NA
        }
        FWD.Cm=FWD.Chm
        REV.Cm=REV.Chm
        FWD.Chm[t(FWD.I) %in% Chm.ids]=Chm.scores[Chm.ids %in% t(FWD.I)]
        REV.Chm[t(REV.I) %in% Chm.ids]=Chm.scores[Chm.ids %in% t(REV.I)]
        FWD.Cm[t(FWD.I) %in% Cm.ids]=Cm.scores[Cm.ids %in% t(FWD.I)]
        REV.Cm[t(REV.I) %in% Cm.ids]=Cm.scores[Cm.ids %in% t(REV.I)]
        FWD.Chm=t(FWD.Chm)
        REV.Chm=t(REV.Chm)
        FWD.Cm=t(FWD.Cm)
        REV.Cm=t(REV.Cm)
        if(length(REV.match)>0 & length(c(REV.match,FWD.match))>1){
          FWD.Chm[(length(FWD.match)+1):(length(REV.match)+length(FWD.match)),]=NA
          FWD.Cm[(length(FWD.match)+1):(length(REV.match)+length(FWD.match)),]=NA
        }
        if(length(FWD.match)>0 & length(c(REV.match,FWD.match))>1){
          REV.Chm[1:length(FWD.match),]=NA
          REV.Cm[1:length(FWD.match),]=NA
        }
        if(length(REV.match)>0 & length(c(REV.match,FWD.match))==1){
          FWD.Chm[]=NA
          FWD.Cm[]=NA
        }
        if(length(FWD.match)>0 & length(c(REV.match,FWD.match))==1){
          REV.Chm[]=NA
          REV.Cm[]=NA
        }
        colnames(FWD.Chm)<-paste0(FWD,'.',1:length(FWD))
        colnames(REV.Chm)<-paste0(REV,'.',1:length(REV))
        colnames(FWD.Cm)<-paste0(FWD,'.',1:length(FWD))
        colnames(REV.Cm)<-paste0(REV,'.',1:length(REV))
        #
        #
        DATA=list('FWD.align'=FWD.align,'REV.align'=REV.align,'FWD.Chm'=FWD.Chm,'FWD.Cm'=FWD.Cm,'REV.Chm'=REV.Chm,'REV.Cm'=REV.Cm,'Q'=quality,'FWD.I'=FWD.I,'REV.I'=REV.I,'Nfrag'=length(c(FWD.match,REV.match)))
      }else{
        DATA=list('Nfrag'=0)
      }
    }
  }
  return(DATA)
}

# Kmer mapping algorithm
map.kmers <- function(read,Cm,Chm,C.key,read.length,FWD,REV,fwdK,revK,Fkey,Rkey) {
  
  k.merge <- function(data){
    
    find.indel <- function(data.sort){
      res=rep(F,times=nrow(data.sort))
      for(i in 1:(nrow(data.sort)-1)){
        if(sum(which(!is.na(data.sort[i,])) %in% which(!is.na(data.sort[i+1,])))==0){
          if(mean(which(!is.na(c(data.sort[i,])))) < mean(which(!is.na(c(data.sort[i+1,]))))){
            if(diff(range(na.omit(c(data.sort[i:(i+1),]))))<=(ncol(data.sort)+ins.tol)){
              res[i]=T
            }
          }
        }
      }
      return(res)
    }
    
    data.I=as.matrix(data[order(rowMeans(data,na.rm = T)),])
    if(length(which(rowMeans(abs(diff(data.I)))==0))>0){
      data.I=as.matrix(data.I[-which(rowMeans(abs(diff(data.I)))==0),])
      if(ncol(data.I)<length(FWD)){
        data.I=t(data.I)
      }
      if(nrow(data.I)==1){
        return(data.I)
      }
    }
    #
    indel.id=which(find.indel(data.I))
    #
    if(length(indel.id)==0){
      return(data.I)
    }else{
      for(i in length(indel.id):1){
        data.I[indel.id[i],which(!is.na(data.I[indel.id[i]+1,]))]=data.I[indel.id[i]+1,which(!is.na(data.I[indel.id[i]+1,]))]
      }
      data.If=data.I[-(indel.id+1),]
      #
      return(data.If)
    }
  }
  
  k.bounds <- function(data){
    slideSum <- function(x,n=1){
      data=c(rep(0,times=n),x,rep(0,times=n))
      results=rep(NA,times=length(x))
      for(i in 1:length(x)){
        results[i]=sum(data[i:(i+2*n)])
      }
      return(results)
    }
    s1=slideSum(as.integer(data))
    s2=range(which(s1>1 & as.numeric(data)))
    results=rep(F,times=length(data))
    results[s2[1]:s2[2]]=T
    return(results)
  }
  
  if(nchar(read)<min(read.length) | nchar(read)>max(read.length)){
    DATA='blank' # bypasses polymers outside desired length range
  } else {
    if((length(Cm)+length(Chm))!=length(C.key)){
      warning('WARNING! Read MM and ML tag info are different lengths...skipping read!')
      DATA='blank'
    }else{
      ImapF=dplyr::bind_rows(lapply(FUN = regexp,X = fwdK,s = read))
      if(length(ImapF)>0){
        if(dim(ImapF)[1]==1){
          Frun = matrix(c(ImapF$start,Fkey[ImapF$match]),2,1)
          Froots = 1
        }else{
          Frun <- rbind(ImapF$start,Fkey[ImapF$match])[,order(ImapF$start)]
          Fbreak = which(diff((Frun[1,]))>=(length(FWD)) | diff((Frun[2,]))<0)
          Findel = which(diff(Frun[1,])!=diff(Frun[2,]) & abs(diff(Frun[1,])-diff(Frun[2,]))<3)
          Findel = Findel[!(Findel %in% Fbreak)]
          Froots = sort(unique(c(Findel,Fbreak,ncol(Frun))))
        }
      }else{
        Findel=rep(0,times=0)
        Froots=rep(0,times=0)
      }
      ImapR=dplyr::bind_rows(lapply(FUN = regexp,X = revK,s = read))
      if(length(ImapR)>0){
        if(dim(ImapR)[1]==1){
          Rrun = matrix(c(ImapR$start,Rkey[ImapR$match]),2,1)
          Rroots = 1
        }else{
          Rrun <- rbind(ImapR$start,Rkey[ImapR$match])[,order(ImapR$start)]
          Rbreak = which(diff(Rrun[1,])>=(length(REV)) | diff(Rrun[2,])<0)
          Rindel = which(diff(Rrun[1,])!=diff(Rrun[2,]) & abs(diff(Rrun[1,])-diff(Rrun[2,]))<3)
          Rindel = Rindel[!(Rindel %in% Rbreak)]
          Rroots = sort(unique(c(Rindel,Rbreak,ncol(Rrun))))
        }
      }else{
        Rindel=rep(0,times=0)
        Rroots=rep(0,times=0)
      }
      if((length(Froots)+length(Rroots))==0){
        return(list('Nfrag'=0))
      }
      #
      polyK=str_split(read,'')[[1]]
      #
      if((length(Froots))==0){
        FWD.I=rep(0,times=0)
        nF=0
      }else{
        FWD.I = t(mapply(seq,Frun[1,Froots]-Frun[2,Froots]+1,Frun[1,Froots]-Frun[2,Froots]+length(FWD),SIMPLIFY = T))
        FWD.I[FWD.I<1 | FWD.I > length(polyK)]=NA
        FWD.I[!t(apply(X = (matrix(polyK[FWD.I],nrow = nrow(FWD.I),ncol = ncol(FWD.I))==matrix(FWD,nrow = nrow(FWD.I),ncol = ncol(FWD.I),byrow = T)),MARGIN = 1,FUN = k.bounds))]=NA
        if(nrow(FWD.I)>1){
          FWD.I=as.matrix(k.merge(FWD.I))
          if(ncol(FWD.I)<length(FWD)){
            FWD.I=t(FWD.I)
          }
        }
        nF=nrow(FWD.I)
      }
      if((length(Rroots))==0){
        REV.I=rep(0,times=0)
        nR=0
      }else{
        REV.I = t(mapply(seq,Rrun[1,Rroots]-Rrun[2,Rroots]+1,Rrun[1,Rroots]-Rrun[2,Rroots]+length(REV),SIMPLIFY = T))
        REV.I[REV.I<1 | REV.I > length(polyK)]=NA
        REV.I[!t(apply(X = (matrix(polyK[REV.I],nrow = nrow(REV.I),ncol = ncol(REV.I))==matrix(REV,nrow = nrow(REV.I),ncol = ncol(REV.I),byrow = T)),MARGIN = 1,FUN = k.bounds))]=NA
        if(nrow(REV.I)>1){
          REV.I=as.matrix(k.merge(REV.I))
          if(ncol(REV.I)<length(REV)){
            REV.I=t(REV.I)
          }
        }
        nR=nrow(REV.I)
      }
      #
      Chm.ids=which(polyK=='C')[Chm]
      Cm.ids=which(polyK=='C')[Cm]
      Chm.scores=C.key[1:length(Chm)]
      Cm.scores=C.key[(length(Chm)+1):(length(Chm)+length(Cm))]
      #
      if(nR==0){
        REV.align=matrix('x',nrow=nF,ncol = ncol(FWD.I))
        REV.I=matrix(NA,nrow=nF,ncol = ncol(FWD.I))
        REV.Chm=REV.I
        REV.Cm=REV.I
        #
        FWD.align=matrix(polyK[FWD.I],nrow = nF,ncol = ncol(FWD.I)); FWD.align[is.na(FWD.align)]='x'
        FWD.Chm=REV.Chm
        if(modk.fillC){
          FWD.Chm[,FWD=='C']=0; FWD.Chm[is.na(FWD.I)]=NA
        }
        FWD.Cm=FWD.Chm
        FWD.Chm[FWD.I %in% Chm.ids]=Chm.scores[Chm.ids %in% FWD.I]
        FWD.Cm[FWD.I %in% Cm.ids]=Cm.scores[Cm.ids %in% FWD.I]
        #
        quality=cbind(rowSums(FWD.align!='x')/ncol(FWD.align),rowSums(FWD.align==matrix(FWD,nrow(FWD.align),ncol(FWD.align),T))/rowSums(FWD.align!='x'),rep('FWD',times=nF+nR))
      }
      if(nF==0){
        FWD.align=matrix('x',nrow=nR,ncol = ncol(REV.I))
        FWD.I=matrix(NA,nrow=nR,ncol = ncol(REV.I))
        FWD.Chm=FWD.I
        FWD.Cm=FWD.I
        #
        REV.align=matrix(polyK[REV.I],nrow = nR,ncol = ncol(REV.I)); REV.align[is.na(REV.align)]='x'
        REV.Chm=FWD.Chm
        if(modk.fillC){
          REV.Chm[,REV=='C']=0; REV.Chm[is.na(REV.I)]=NA
        }
        REV.Cm=REV.Chm
        REV.Chm[REV.I %in% Chm.ids]=Chm.scores[Chm.ids %in% REV.I]
        REV.Cm[REV.I %in% Cm.ids]=Cm.scores[Cm.ids %in% REV.I]
        #
        quality=cbind(rowSums(REV.align!='x')/ncol(REV.align),rowSums(REV.align==matrix(REV,nrow(REV.align),ncol(REV.align),T))/rowSums(REV.align!='x'),rep('REV',times=nF+nR))
      }
      if(nF>0 & nR>0){
        FWD.I=rbind(FWD.I,matrix(NA,nrow=nR,ncol = ncol(REV.I)))
        REV.I=rbind(matrix(NA,nrow=nF,ncol = ncol(FWD.I)),REV.I)
        #
        FWD.align=matrix(polyK[FWD.I],nrow = nF+nR,ncol = ncol(FWD.I)); FWD.align[is.na(FWD.align)]='x'
        REV.align=matrix(polyK[REV.I],nrow = nF+nR,ncol = ncol(REV.I)); REV.align[is.na(REV.align)]='x'
        FWD.Chm=matrix(NA,nrow=nF+nR,ncol = ncol(REV.I))
        REV.Chm=FWD.Chm
        if(modk.fillC){
          REV.Chm[,REV=='C']=0; REV.Chm[is.na(REV.I)]=NA
          FWD.Chm[,FWD=='C']=0; FWD.Chm[is.na(FWD.I)]=NA
        }
        FWD.Cm=FWD.Chm
        REV.Cm=REV.Chm
        FWD.Chm[FWD.I %in% Chm.ids]=Chm.scores[Chm.ids %in% FWD.I]
        FWD.Cm[FWD.I %in% Cm.ids]=Cm.scores[Cm.ids %in% FWD.I]
        REV.Chm[REV.I %in% Chm.ids]=Chm.scores[Chm.ids %in% REV.I]
        REV.Cm[REV.I %in% Cm.ids]=Cm.scores[Cm.ids %in% REV.I]
        #
        comp=rowSums(cbind(FWD.align,REV.align)!='x')/length(FWD)
        match=rowSums(cbind(FWD.align,REV.align)==matrix(c(FWD,REV),nrow(FWD.align),2*length(FWD),T),na.rm = T)/comp/length(FWD)
        quality=cbind(comp,match,rep(c('FWD','REV'),times=c(nF,nR)))
      }
      DATA=list('FWD.align'=FWD.align,'REV.align'=REV.align,'FWD.Chm'=FWD.Chm,'FWD.Cm'=FWD.Cm,'REV.Chm'=REV.Chm,'REV.Cm'=REV.Cm,'Q'=quality,'FWD.I'=FWD.I,'REV.I'=REV.I,'Nfrag'=nF+nR)
    }
  }
  return(DATA)
}

# Fragment mapping algorithm
map.fragments <- function(read,Cm,Chm,C.key,read.length,FWD,REV) {

  bounds <- function(data){
    slideSum <- function(x,n=1){
      data=c(rep(0,times=n),x,rep(0,times=n))
      results=rep(NA,times=length(x))
      for(i in 1:length(x)){
        results[i]=sum(data[i:(i+2*n)])
      }
      return(results)
    }
    s1=slideSum(as.integer(data))
    if(sum(s1>1 & as.numeric(data))==0){
      return(NULL)
    }
    s2=range(which(s1>1 & as.numeric(data)))
    results=rep(F,times=length(data))
    results[s2[1]:s2[2]]=T
    return(results)
  }

  get.scores <- function(poly.ref,FWD,REV){
    FWD.score=rep(0,times=length(poly.ref)+length(FWD))
    REV.score=rep(0,times=length(poly.ref)+length(REV))
    for(i in 1:length(REV)){
      FWD.score[which(poly.ref==FWD[i])-i+1+length(FWD)]=FWD.score[which(poly.ref==FWD[i])-i+1+length(FWD)]+(1/length(FWD))
      FWD.score[which(poly.ref!=FWD[i] & poly.ref!='x')-i+1+length(FWD)]=FWD.score[which(poly.ref!=FWD[i] & poly.ref!='x')-i+1+length(FWD)]-(1/length(FWD)/10)
      REV.score[which(poly.ref==REV[i])-i+1+length(REV)]=REV.score[which(poly.ref==REV[i])-i+1+length(REV)]+(1/length(REV))
      REV.score[which(poly.ref!=REV[i] & poly.ref!='x')-i+1+length(REV)]=REV.score[which(poly.ref!=REV[i] & poly.ref!='x')-i+1+length(REV)]-(1/length(REV)/10)
    }
    res=data.frame('fwd'=FWD.score,'rev'=REV.score)
    return(res)
  }

  if(nchar(read)<min(read.length) | nchar(read)>max(read.length)){
    DATA='blank' # bypasses polymers outside desired length range
  } else {
    if((length(Cm)+length(Chm))!=length(C.key)){
      warning('WARNING! Read MM and ML tag info are different lengths...skipping read!')
      DATA='blank'
    }else{
      #
      NN=round(nchar(read)/mean(c(length(FWD),length(REV)))+0.4)+padding # sets the maximum number of fragments allowed to be mapped to the polymer
      poly.ref=str_split(read,'')[[1]]
      Chm.ids=which(poly.ref=='C')[Chm]
      Cm.ids=which(poly.ref=='C')[Cm]
      Chm.scores=C.key[1:length(Chm)]
      Cm.scores=C.key[(length(Chm)+1):(length(Chm)+length(Cm))]
      COUNTER=1
      #
      FWD.align=matrix('x',nrow=NN,ncol = length(FWD)); colnames(FWD.align)<-paste0(FWD,'.',1:length(FWD))
      REV.align=matrix('x',nrow=NN,ncol = length(REV)); colnames(REV.align)<-paste0(REV,'.',1:length(REV))
      FWD.Chm=matrix(NA,nrow=NN,ncol = length(FWD)); colnames(FWD.Chm)<-paste0(FWD,'.',1:length(FWD))
      REV.Chm=matrix(NA,nrow=NN,ncol = length(REV)); colnames(REV.Chm)<-paste0(REV,'.',1:length(REV))
      FWD.Cm=FWD.Chm
      REV.Cm=REV.Chm
      FWD.I=matrix(NA,nrow=NN,ncol = length(FWD))
      REV.I=matrix(NA,nrow=NN,ncol = length(REV))
      Q=as.data.frame(matrix(NA,nrow=NN,ncol = 3));colnames(Q) <- c('comp','match','id')
      #
      perf.matches=find.motif(read,Cm,Chm,C.key,read.length,FWD,REV)
      if(perf.matches[['Nfrag']]==1){
        FWD.align[COUNTER,]=perf.matches[['FWD.align']]
        REV.align[COUNTER,]=perf.matches[['REV.align']]
        FWD.Chm[COUNTER,]=perf.matches[['FWD.Chm']]
        REV.Chm[COUNTER,]=perf.matches[['REV.Chm']]
        FWD.Cm[COUNTER,]=perf.matches[['FWD.Cm']]
        REV.Cm[COUNTER,]=perf.matches[['REV.Cm']]
        FWD.I[COUNTER,]=perf.matches[['FWD.I']]
        REV.I[COUNTER,]=perf.matches[['REV.I']]
        Q[COUNTER,]=perf.matches[['Q']]
        #
        poly.ref[as.numeric(c(perf.matches[['FWD.I']][perf.matches[['FWD.I']]>0],perf.matches[['REV.I']][perf.matches[['REV.I']]>0]))]='x'
        COUNTER=COUNTER+1
      }
      if(perf.matches[['Nfrag']]>1){
        FWD.align[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['FWD.align']]
        REV.align[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['REV.align']]
        FWD.Chm[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['FWD.Chm']]
        REV.Chm[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['REV.Chm']]
        FWD.Cm[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['FWD.Cm']]
        REV.Cm[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['REV.Cm']]
        FWD.I[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['FWD.I']]
        REV.I[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['REV.I']]
        Q[COUNTER:(COUNTER+perf.matches[['Nfrag']]-1),]=perf.matches[['Q']]
        #
        poly.ref[as.numeric(c(perf.matches[['FWD.I']][perf.matches[['FWD.I']]>0],perf.matches[['REV.I']][perf.matches[['REV.I']]>0]))]='x'
        COUNTER=COUNTER+perf.matches[['Nfrag']]
      }
      #
      if(COUNTER<=NN){
        for (j in COUNTER:NN){
          if(sum(poly.ref!="x")<10){
            break # terminates further alignment attempts when less than 10 bases in the polymer remain unmapped
          }
          #
          align.scores=get.scores(poly.ref,FWD,REV)
          if(max(align.scores$fwd)>=max(align.scores$rev)){
            edges=bounds(c(rep('x',times=length(FWD)),poly.ref,rep('x',times=length(FWD)))[(which.max(align.scores$fwd)):(which.max(align.scores$fwd)+length(FWD)-1)]==FWD)
            if(is.null(edges)){
              break
            }
            FWD.align[j,which(edges)]=poly.ref[(which.max(align.scores$fwd)-length(FWD)+which(edges)-1)]
            if(modk.fillC){
              FWD.Chm[j,which(FWD.align[j,]=='C')]=0; FWD.Chm[j,which(!edges)]=NA
              FWD.Cm[j,which(FWD.align[j,]=='C')]=0; FWD.Cm[j,which(!edges)]=NA
            }
            FWD.Chm[j,which(edges)[(which.max(align.scores$fwd)-length(FWD)+which(edges)-1) %in% Chm.ids]]=Chm.scores[Chm.ids %in% (which.max(align.scores$fwd)-length(FWD)+which(edges)-1)]
            FWD.Cm[j,which(edges)[(which.max(align.scores$fwd)-length(FWD)+which(edges)-1) %in% Cm.ids]]=Cm.scores[Cm.ids %in% (which.max(align.scores$fwd)-length(FWD)+which(edges)-1)]
            FWD.I[j,which(edges)]=(which.max(align.scores$fwd)-length(FWD)+which(edges)-1)
            #
            Q[j,1]=sum(edges)/length(FWD)
            Q[j,2]=mean(poly.ref[(which.max(align.scores$fwd)-length(FWD)+which(edges)-1)]==FWD[edges])
            Q[j,3]='FWD'
            #
            poly.ref[(which.max(align.scores$fwd)-length(FWD)+which(edges)-1)]='x'
          }
          if(max(align.scores$fwd)<max(align.scores$rev)){
            edges=bounds(c(rep('x',times=length(REV)),poly.ref,rep('x',times=length(REV)))[(which.max(align.scores$rev)):(which.max(align.scores$rev)+length(REV)-1)]==REV)
            if(is.null(edges)){
              break
            }
            REV.align[j,which(edges)]=poly.ref[(which.max(align.scores$rev)-length(REV)+which(edges)-1)]
            if(modk.fillC){
              REV.Chm[j,which(REV.align[j,]=='C')]=0; REV.Chm[j,which(!edges)]=NA
              REV.Cm[j,which(REV.align[j,]=='C')]=0; REV.Cm[j,which(!edges)]=NA
            }
            REV.Chm[j,which(edges)[((which.max(align.scores$rev)-length(REV)+which(edges)-1) %in% Chm.ids)]]=Chm.scores[Chm.ids %in% (which.max(align.scores$rev)-length(REV)+which(edges)-1)]
            REV.Cm[j,which(edges)[((which.max(align.scores$rev)-length(REV)+which(edges)-1) %in% Cm.ids)]]=Cm.scores[Cm.ids %in% (which.max(align.scores$rev)-length(REV)+which(edges)-1)]
            REV.I[j,which(edges)]=(which.max(align.scores$rev)-length(REV)+which(edges)-1)
            #
            Q[j,1]=sum(edges)/length(REV)
            Q[j,2]=mean(poly.ref[(which.max(align.scores$rev)-length(REV)+which(edges)-1)]==REV[edges])
            Q[j,3]='REV'
            #
            poly.ref[(which.max(align.scores$rev)-length(REV)+which(edges)-1)]='x'
          }
        }
      }

      # save relevant data from read to subset of stable variable
      DATA=list('FWD.align'=FWD.align,'REV.align'=REV.align,'FWD.Chm'=FWD.Chm,'FWD.Cm'=FWD.Cm,'REV.Chm'=REV.Chm,'REV.Cm'=REV.Cm,'Q'=as.matrix(Q),'FWD.I'=FWD.I,'REV.I'=REV.I)
      if((sum(DATA$Q[,3]=='FWD',na.rm = T)>1 | sum(DATA$Q[,3]=='REV',na.rm = T)>1)){
        DATA=indel.fix(DATA,FWD,REV)
      }

    }
  }

  return(DATA)

}

# Mixed beta thresholding
BM.thresh <- function(data.actual.fwd,data.actual.rev,met='RSS',p=(1-T1.err),set=BM.strand.data,p2=T2.err){
  if(set=='b'){
    data=as.numeric(na.omit(c(data.actual.fwd[data.actual.fwd>0 & data.actual.fwd<1],data.actual.rev[data.actual.rev>0 & data.actual.rev<1])))
  }
  if(set=='f'){
    data=as.numeric(na.omit(c(data.actual.fwd[data.actual.fwd>0 & data.actual.fwd<1])))
  }
  if(set=='r'){
    data=as.numeric(na.omit(c(data.actual.rev[data.actual.rev>0 & data.actual.rev<1])))
  }
  #i=1;data=as.numeric(na.omit(c(data.actual.rev[data.actual.rev[,i]>0 & data.actual.rev[,i]<1,i])))
  fit.dens=density(x = data,na.rm = T,from = 0,to = 1,width = 0.05);fit.dens$y=fit.dens$y/sum(fit.dens$y)/mean(diff(fit.dens$x))
  beta.par.calc <- function(m,s){
    alpha=m*(m*(1-m)/s^2-1)
    beta=alpha*(1-m)/m
    return(as.numeric(c(alpha,beta)))
  }
  reg.beta <- function(par,dens=fit.dens,dat=data,meth=met){
    if(meth=='NLL'){
      NLL=-sum(log(par[5]*dbeta(dat,shape1 = par[1],shape2 = par[2])+(1-par[5])*dbeta(dat,shape1 = par[3],shape2 = par[4])))
      return(NLL)
    }
    if(meth=='RSS'){
      b.pdf=par[5]*dbeta(dens$x[-c(1,length(dens$x))],shape1 = par[1],shape2 = par[2])+(1-par[5])*dbeta(dens$x[-c(1,length(dens$x))],shape1 = par[3],shape2 = par[4])
      res=sum((dens$y[-c(1,length(dens$x))]-b.pdf)^2)
      return(res)
    }
  }
  init.par.est <- function(data,fit.dens){
    sep.point=fit.dens$x[round(mean(pracma::findpeaks(fit.dens$y,minpeakheight = max(fit.dens$y,na.rm = T)*0.2,sortstr = T,npeaks = 2)[,2]))]
    p=sum(fit.dens$y[fit.dens$x<sep.point],na.rm = T)/sum(fit.dens$y,na.rm = T)
    m1=mean(data[data<sep.point],na.rm = T)
    sd1=sd(data[data<sep.point],na.rm = T)
    m2=mean(data[data>=sep.point],na.rm = T)
    sd2=sd(data[data>=sep.point],na.rm = T)
    pars=c(a1=beta.par.calc(m1,sd1)[1],b1=beta.par.calc(m1,sd1)[2],a2=beta.par.calc(m2,sd2)[1],b2=beta.par.calc(m2,sd2)[2],p=p)
    return(pars)
  }
  reg.beta.single <- function(par,dens=fit.dens,dat=data,meth=met){
    if(meth=='NLL'){
      NLL=-1*sum(log(dbeta(dat,shape1 = par[1],shape2 = par[2])))
      return(NLL)
    }
    if(meth=='RSS'){
      b.pdf=dbeta(dens$x[-c(1,length(dens$x))],shape1 = par[1],shape2 = par[2])
      res=sum((dens$y[-c(1,length(dens$x))]-b.pdf)^2)
      return(res)
    }
  }
  init.par.est.single <- function(data,fit.dens){
    m=mean(data)
    sd=sd(data)
    pars=c(a=beta.par.calc(m,sd)[1],b=beta.par.calc(m,sd)[2])
    return(pars)
  }
  try(fit.betasSTD <- optim(par = c(a1=7,b1=1.2,a2=5,b2=50,p=0.5),fn = reg.beta,lower = c(0,0,0,0,0),upper = c(Inf,Inf,Inf,Inf,1),method = "L-BFGS-B"))
  try(fit.betasEST <- optim(par = init.par.est(data,fit.dens),fn = reg.beta,lower = c(0,0,0,0,0),upper = c(Inf,Inf,Inf,Inf,1),method = "L-BFGS-B"))
  try(fit.betasSIN <- optim(par = init.par.est.single(data,fit.dens),fn = reg.beta.single))
  if(exists('fit.betasEST') & !exists('fit.betasSTD')){
    fit.betas=fit.betasEST
  }
  if(!exists('fit.betasEST') & exists('fit.betasSTD')){
    fit.betas=fit.betasSTD
  }
  if(exists('fit.betasEST') & exists('fit.betasSTD')){
    if(fit.betasEST$value<=fit.betasSTD$value){
      fit.betas=fit.betasEST
    }
    if(fit.betasEST$value>fit.betasSTD$value){
      fit.betas=fit.betasSTD
    }
  }
  if(!exists('fit.betasSIN')){
    try(fit.betasSINn <- optim(par = c(a=8,b=100),fn = reg.beta.single))
    try(fit.betasSINp <- optim(par = c(a=100,b=8),fn = reg.beta.single))
    if(fit.betasSINn$value < fit.betasSINp$value){
      fit.betasSIN=fit.betasSINn
    }
    if(fit.betasSINp$value <= fit.betasSINn$value){
      fit.betasSIN=fit.betasSINp
    }
  }
  dBeta.share=sum(apply(X = rbind(fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]),(1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4])),MARGIN = 2,FUN = min))*mean(diff(fit.dens$x),na.rm=T)/min(c(fit.betas$par[5],1-fit.betas$par[5]))
  dBIC=2*length(fit.betas$par)+length(fit.dens$x)*log(fit.betas$value/length(fit.dens$x))
  sBIC=2*length(fit.betasSIN$par)+length(fit.dens$x)*log(fit.betasSIN$value/length(fit.dens$x))
  mode.diff=abs(fit.dens$x[which.max(dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]))]-fit.dens$x[which.max(dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]))])
  if(dBIC<sBIC & dBeta.share<0.5 & abs(0.5-fit.betas$par[5])<0.4 & mode.diff>0.1){
    beta.means=as.numeric((fit.betas$par[c(1,3)])/(fit.betas$par[c(1,3)]+fit.betas$par[c(2,4)]))
    if(which.max(beta.means)==1){
      rel.lik=((fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]))/((1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4])))
    }
    if(which.max(beta.means)==2){
      rel.lik=1/((fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]))/((1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4])))
    }
    thresh=qbeta(p,fit.betas$par[2*(which.min(beta.means)-1)+1],fit.betas$par[2*(which.min(beta.means)-1)+2])
    #thresh2=fit.dens$x[min(which(rel.lik>1 & fit.dens$x>min(beta.means)))]
    plot(fit.dens,ylim=c(0,max(fit.dens$y)),main='Beta Unmixing Threshold',xlab='Methyl Score')
    lines(fit.dens$x,fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2])+(1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]),col='red')
    abline(v=thresh,col='green',lwd=2,lty='dashed')
    #abline(v=thresh2,col='blue',lwd=2,lty='dashed')
    legend('topright',legend = c('Data','Model','Threshold'),col=c('black','red','green'),fill=c('black','red','green'))
    lines(fit.dens$x,fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]),lty='dotted',col='purple',lwd=2)
    lines(fit.dens$x,(1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]),lty='dotted',col='purple',lwd=2)
    text(x=1,y=0.85*max(fit.dens$y),pos=2,labels=paste0('T1 ≈ ',round(1-p,2),' (',round(100*(1-p)*c(fit.betas$par[5],1-fit.betas$par[5])[which.min(beta.means)]),'%)'),cex=2,col='black')
    false.neg=pbeta(thresh,fit.betas$par[2*(which.max(beta.means)-1)+1],fit.betas$par[2*(which.max(beta.means)-1)+2])
    text(x=1,y=0.75*max(fit.dens$y),pos=2,labels=paste0('T2 ≈ ',round(false.neg,2),' (',round(100*false.neg*c(fit.betas$par[5],1-fit.betas$par[5])[which.max(beta.means)]),'%)'),cex=2,col='black')
  }
  if(sBIC<dBIC | dBeta.share>=0.5 | abs(0.5-fit.betas$par[5])>0.4 | mode.diff<=0.1){
    plot(fit.dens,ylim=c(0,max(fit.dens$y)),main='Beta Unmixing Threshold',xlab='Methyl Score');lines(fit.dens$x,dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]),col='red')
    show('WARNING! -- Methyl score values appear to primarily belong to a single distribution, and its corresponding methylation state is needed for thresholding -- attempting auto-assignment...')
    if(fit.dens$x[which.max(fit.dens$y)]>=0.3){
      usr.input='p'
      show('...distribution inferred to correspond to methylated CpGs.')
    }else{
      if(fit.dens$x[which.max(fit.dens$y)]<=0.1){
        usr.input='n'
        show('...distribution inferred to correspond to unmethylated CpGs.')
      }else{
        dist.diff=fit.dens$y-dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]);dist.diff[dist.diff<0]=0
        if(sum(dist.diff[1:which.max(fit.dens$y)],na.rm = T)<sum(dist.diff[which.max(fit.dens$y):length(fit.dens$y)],na.rm = T)){
          usr.input='n'
          show('...distribution inferred to correspond to unmethylated CpGs.')
        }
        if(sum(dist.diff[1:which.max(fit.dens$y)],na.rm = T)>=sum(dist.diff[which.max(fit.dens$y):length(fit.dens$y)],na.rm = T)){
          usr.input='p'
          show('...distribution inferred to correspond to methylated CpGs.')
        }
      }
    }
    if(usr.input=='n'){
      thresh=qbeta(p,fit.betasSIN$par[1],fit.betasSIN$par[2])
      pos.dens=fit.dens$y-dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]);pos.dens[pos.dens<0 | fit.dens$x<fit.dens$x[which.max(dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]))]]=0
      if(sum(pos.dens*mean(diff(fit.dens$x),na.rm=T))>=0.05){
        false.neg=(cumsum(pos.dens)/sum(pos.dens))[which.min(abs(fit.dens$x-thresh))]
        text(x=1,y=0.85*max(fit.dens$y),pos=2,labels=paste0('T1 ≈ ',round(1-p,2),' (',round(100*(1-p)*(1-sum(pos.dens*mean(diff(fit.dens$x),na.rm=T)))),'%)'),cex=2,col='black')
        text(x=1,y=0.75*max(fit.dens$y),pos=2,labels=paste0('T2 ≈ ',round(false.neg,2),' (',round(100*false.neg*sum(pos.dens*mean(diff(fit.dens$x),na.rm=T))),'%)'),cex=2,col='black')
        lines(fit.dens$x,pos.dens,lty='dotted',col='purple',lwd=2)
        lines(fit.dens$x,(1-mean(diff(fit.dens$x),na.rm=T)*sum(pos.dens,na.rm = T))*dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]),lty='dotted',col='purple',lwd=2)
      }else{
        text(x=1,y=0.85*max(fit.dens$y),pos=2,labels=paste0('T1 ≈ ',round(1-p,2),' (≤',round(100*(1-p)),'%)'),cex=2,col='black')
        text(x=1,y=0.75*max(fit.dens$y),pos=2,labels=paste0('T2 ≈ ','n.d.',' (','≤5','%)'),cex=2,col='black')
      }
    }
    if(usr.input=='p'){
      neg.dens=fit.dens$y-dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]);neg.dens[neg.dens<0 | fit.dens$x<fit.dens$x[which.max(dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]))]]=0
      if(sum(neg.dens*mean(diff(fit.dens$x),na.rm=T))>=0.05){
        thresh=fit.dens$x[min(which((cumsum(neg.dens)/sum(neg.dens))>p))]
        text(x=1,y=0.85*max(fit.dens$y),pos=2,labels=paste0('T1 ≈ ',round(1-p,2),' (',round(100*(1-p)*sum(neg.dens*mean(diff(fit.dens$x),na.rm=T))),'%)'),cex=2,col='black')
        text(x=1,y=0.75*max(fit.dens$y),pos=2,labels=paste0('T2 ≈ ',round(pbeta(thresh,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]),2),' (',round(100*pbeta(thresh,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2])*(1-sum(neg.dens*mean(diff(fit.dens$x),na.rm=T)))),'%)'),cex=2,col='black')
        lines(fit.dens$x,neg.dens,lty='dotted',col='purple',lwd=2)
        lines(fit.dens$x,(1-mean(diff(fit.dens$x),na.rm=T)*sum(neg.dens,na.rm = T))*dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]),lty='dotted',col='purple',lwd=2)
      }else{
        thresh=qbeta(p2,fit.betasSIN$par[1],fit.betasSIN$par[2])
        text(x=1,y=0.85*max(fit.dens$y),pos=2,labels=paste0('T1 ≈ ','n.d.',' (','≤5','%)'),cex=2,col='black')
        text(x=1,y=0.75*max(fit.dens$y),pos=2,labels=paste0('T2 ≈ ',round(p2,2),' (≤',round(100*p2),'%)'),cex=2,col='black')
      }
    }
    abline(v=thresh,col='green',lwd=2,lty='dashed')
    legend('topright',legend = c('Data','Model','Threshold'),col=c('black','red','green'),fill=c('black','red','green'))
  }
  return(thresh)
}

# scanning function to determine binarization threshold value
find_thresh_for_auc <- function(data.actual.fwd, data.actual,target_auc = 0.9, n_initial = 30, refine_steps = 5) {

  # helper function for heterogenous read length analysis
  trim_trailing_nas <- function(mat) {
    # Trim trailing NAs row-wise
    trimmed_rows <- lapply(seq_len(nrow(mat)), function(i) {
      row <- mat[i, ]
      if (all(is.na(row))) return(NA_real_)
      last_non_na <- max(which(!is.na(row)))
      row[1:last_non_na]
    })
    # Determine max row length after trimming
    max_len <- max(lengths(trimmed_rows))
    # Pad rows back to rectangular matrix
    padded <- t(sapply(trimmed_rows, function(row) {
      c(row, rep(NA_real_, max_len - length(row)))
    }))
    storage.mode(padded) <- "numeric"
    return(padded)
  }

  #trim trailing NAs & repad to matrix
  fwd_trimmed <- trim_trailing_nas(data.actual.fwd)
  rev_trimmed <- trim_trailing_nas(data.actual)
  # flatten numeric range for threshold search
  flat <- unlist(fwd_trimmed)
  flat <- as.numeric(flat[!is.na(suppressWarnings(as.numeric(flat)))])
  search_range <- range(flat, na.rm = TRUE)
  # AUC calculator for a given threshold
  compute_fwd_auc <- function(thresh) {
    # binarize  matrices
    fwd.binary <- binarize_matrix(fwd_trimmed, thresh)
    rev.binary <- binarize_matrix(rev_trimmed, thresh)
    # reverse columns (cat strand)
    rev.binary <- rev.binary[, rev(seq_len(ncol(rev.binary))), drop = FALSE]
    # compute column means
    fwd.binary.mean <- colMeans(fwd.binary, na.rm = TRUE)
    # keep positions with appreciable forward methylation
    keep <- which(fwd.binary.mean > 0.5)
    if (length(keep) < 2) return(NA)
    fwd.filtered <- fwd.binary[, keep, drop = FALSE]
    # survival durations
    surv_vals <- get_survival_data(fwd.filtered)$duration
    if (length(surv_vals) < 2) return(NA)
    surv_curve <- empirical_survival(surv_vals)
    compute_auc(surv_curve)
  }
  # initial scan for threshold value
  thresh_grid <- seq(search_range[1], search_range[2], length.out = n_initial)
  auc_vals <- sapply(thresh_grid, compute_fwd_auc)
  # refinement
  for (i in seq_len(refine_steps)) {
    valid <- which(!is.na(auc_vals))
    if (length(valid) < 3) break
    idx_best <- valid[which.min(abs(auc_vals[valid] - target_auc))]
    best_thresh <- thresh_grid[idx_best]
    window <- diff(search_range) / (5 * i)
    sub_range <- c(
      max(search_range[1], best_thresh - window),
      min(search_range[2], best_thresh + window)
    )
    thresh_grid <- seq(sub_range[1], sub_range[2], length.out = n_initial)
    auc_vals <- sapply(thresh_grid, compute_fwd_auc)
  }
  #final threshold value selection
  valid <- which(!is.na(auc_vals))
  idx_best <- valid[which.min(abs(auc_vals[valid] - target_auc))]
  best_thresh <- thresh_grid[idx_best]
  final_auc  <- auc_vals[idx_best]
  message(sprintf("Chosen threshold = %.4f (achieved AUC = %.3f)",
                  best_thresh, final_auc))
  return(best_thresh)
}

# binarize data based on threshold (methylated=1, unmethylated=0)
binarize_matrix <- function(mat, threshold) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input must be a matrix or data frame.")
  }
  # convert to numeric matrix
  mat_num <- as.matrix(mat)
  # apply threshold: binarize methylation data
  if(length(threshold)>1){
    bin_mat <- ifelse(mat_num < matrix(threshold,nrow = nrow(mat_num),ncol = ncol(mat_num),byrow = T), 0, 1)
  }else{
    bin_mat <- ifelse(mat_num < threshold, 0, 1)
  }
  return(bin_mat)
}

# survival analysis - find survival processivity metrics
get_survival_data <- function(mat) {
  do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
    row <- mat[i, ]
    # skip if entire row is NA
    if (all(is.na(row))) return(NULL)
    # trim trailing NAs (for heterogenous length products)
    non_na_idx <- which(!is.na(row))
    if (length(non_na_idx) == 0) return(NULL)
    last_non_na <- max(non_na_idx)
    row <- row[1:last_non_na]
    # find the first methylation event
    first_one_idx <- which(row == 1)
    if (length(first_one_idx) == 0) return(NULL)
    first_one <- first_one_idx[1]
    # slice from first methylation onward
    trimmed <- row[first_one:length(row)]
    # rle on trimmed region
    run <- rle(trimmed)
    consecutive_ones <- run$lengths[1]
    value_first <- run$values[1]
    if (value_first != 1) {
      return(data.frame(read_id = i, duration = 0))
    }
    total_len <- length(trimmed)
    duration <- consecutive_ones / total_len
    # single isolated methylation duration always = 0
    if (consecutive_ones == 1) duration <- 0
    data.frame(read_id = i, duration = duration)
  }))
}

# put survival data into appropriate format for survival analysis
empirical_survival <- function(durations) {
  durations <- sort(durations)
  n <- length(durations)
  # count zeros
  zero_count <- sum(durations == 0)
  nonzero_durations <- durations[durations > 0]
  # build step_x and step_y for survival
  step_x <- c(0, 0, nonzero_durations, 1)
  step_y <- c(100, (n - zero_count)/n * 100, (rev(seq_along(nonzero_durations)))/n * 100, 0)
  data.frame(x = step_x, y = step_y)
}

# fraction of reads fully methylated
fraction_full <- function(durations) mean(durations == 1)

# survival analysis AUC calc
compute_auc <- function(surv_df) {
  x <- surv_df$x
  y <- surv_df$y / 100
  auc <- sum(diff(x) * head(y, -1))
  return(auc)
}

# median survival duration calc
median_duration <- function(surv_df) {
  idx <- which(surv_df$y <= 50)[1]
  return(surv_df$x[idx])
}

# row summary stats for distribution analysis
row_sum_counts <- function(mat) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input must be a matrix or data frame.")
  }
  rs <- rowSums(mat)
  total_rows <- length(rs)
  df <- data.frame(
    sum = 0:max(rs),
    count = sapply(0:max(rs), function(k) sum(rs == k))
  )
  df$percentage <- df$count / total_rows * 100
  return(df)
}

# position scoring for mapped fragments
pos.score <- function(indexes,n,k){
  data=indexes;data[indexes<=0]=NA
  midpoint=apply(data,1,median,na.rm=T)
  score=cbind(midpoint/n,(n/2-abs(n/2-midpoint))/k,midpoint/k,(midpoint-n)/k,(midpoint-n/2)/k)
  return(score)
}

# probability density function maker
pdf.make <- function(data,pars=NULL){
  if(is.null(pars)){
    tmp1=density(data,na.rm = T)
  }else{
    if(length(pars)==2){
      tmp1=density(data,from = pars[1],to = pars[2],na.rm = T)
    }
    if(length(pars)==3){
      tmp1=density(data,from = pars[1],to = pars[2],bw = pars[3],na.rm = T)
    }
  }
  tmp1$y=tmp1$y/sum(mean(diff(tmp1$x))*tmp1$y)
  return(tmp1)
}

# generate modkit summary-like statistics
modk.sum <- function(melo,meca,seqq,thresh,methyl.type,fill.Cs){
  if(!fill.Cs){
    meca.set=as.numeric((str_split(meca,pattern = ',')[[1]])[-1])
    Nc=length(meca.set)/2
    Cm.scores=matrix(NA,nrow=2,ncol=Nc)
    Cm.scores[1,]=meca.set[1:Nc]/256
    Cm.scores[2,]=meca.set[(Nc+1):(2*Nc)]/256
    if(methyl.type=='B'){
      res=mean(colSums(Cm.scores)>=thresh)
    }
    if(methyl.type=='M'){
      res=mean(Cm.scores[2,]>=thresh)
    }
    if(methyl.type=='H'){
      res=mean(Cm.scores[1,]>=thresh)
    }
    return(data.frame('pC'=res,'Nc'=Nc))
  }else{
    meloCH.id=cumsum(as.numeric((str_split(melo[,1],pattern = ',')[[1]])[-1])+1)
    meloCM.id=cumsum(as.numeric((str_split(melo[,1],pattern = ',')[[1]])[-1])+1)
    meca.set=as.numeric((str_split(meca,pattern = ',')[[1]])[-1])
    Nc=sum(str_split(seqq,pattern = '')[[1]]=='C',na.rm = T)
    Cm.scores=matrix(0,nrow=2,ncol=Nc)
    Cm.scores[1,meloCH.id]=meca.set[1:length(meloCH.id)]/256
    Cm.scores[2,meloCM.id]=meca.set[(length(meloCH.id)+1):(length(meloCH.id)+length(meloCM.id))]/256
    if(methyl.type=='B'){
      res=mean(colSums(Cm.scores)>=thresh)
    }
    if(methyl.type=='M'){
      res=mean(Cm.scores[2,]>=thresh)
    }
    if(methyl.type=='H'){
      res=mean(Cm.scores[1,]>=thresh)
    }
    return(data.frame('pC'=res,'Nc'=Nc))
  }
}

# Adapter for converting bisulfite sequencing FASTQ files into MeLoELF-compatible input txt files
fastq.adapter <- function(read){
  eM.Z=diff(c(0,which(as.numeric(gregexpr('C|U',read)[[1]]) %in% as.numeric(gregexpr('U',read)[[1]]))))-1
  res=data.frame('read'=gsub('U','C',read),'melo'=paste0('EM:Z:C+h?',paste0(c(rbind(rep(',',times=length(eM.Z)),eM.Z)),collapse = ''),';C+m?',paste0(c(rbind(rep(',',times=length(eM.Z)),eM.Z)),collapse = ''),';',collapse = ''),'meca'=paste0('EL:B:C',paste0(rep(c(',',0),each=length(eM.Z)),collapse = ''),paste0(rep(c(',',256),each=length(eM.Z)),collapse = ''),collapse = ''))
  return(res)
}

# Save relevant parameters
PAR.list=ls()

#######################
## Read Parsing and Fragment Mapping
#######################


if(crunch.too){

  # automatic search for appropriate input file
  if(sam.file=='auto'){
    sam.file=list.files(path = mdir,pattern = '*.txt')
    if(sum(c(seq.file,meca.file,melo.file) %in% sam.file)==3){
      show('TXT files detected...bypassing pre-processing.')
      sam.file=NULL
    }else{
      show('Pre-processed TXT files matching the relevant parameter names are unavailable...searching for SAM file.')
      sam.file=list.files(path = mdir,pattern = '*.sam')
      if(length(sam.file==1)){
        show("SAM file detected!")
      }else{
        show("Set sam.file='auto', but a single SAM file was not found in working directory...defaulting to search for FASTQ file.")
        sam.file=list.files(path = mdir,pattern = '*.fastq')
        if(length(sam.file)==1){
          show('FASTQ file detected...it will be treated as bisulfite sequencing data, and input adapted accordingly. Adapting (~50k reads/min)...')
          used.fastq=T
          if(!modk.fillC){
            stop("ERROR: Bisulfite sequencing FASTQ files should be analyzed with modk.fillC=T to ensure proper methylation assignment for C basecalls.")
          }
        }else{
          show('No single FASTQ file found in working directory...')
          stop('ERROR! No proper input files found -- terminating run.')
        }
      }
    }
  }
  
  # process fastq file into individual txt files for loading
  if(used.fastq){
    pre.q=fread(sam.file,header = F,sep = '\t')
    pre.q=pre.q[seq(2,nrow(pre.q),4),1]
    registerDoParallel(detectCores())
    ONTsim <- foreach(i = 1:length(pre.q),.combine = 'rbind') %dopar% {
      fastq.adapter(pre.q[i])
    }
    fwrite(x = as.list(ONTsim$read),file = seq.file,sep = '\n',col.names = F,quote = F)
    fwrite(x = as.list(ONTsim$melo),file = melo.file,sep = '\n',col.names = F,quote = F)
    fwrite(x = as.list(ONTsim$meca),file = meca.file,sep = '\n',col.names = F,quote = F)
    show('...file adaptation complete!')
  }
  
  # process relevant sam file information into individual txt files for loading
  if(!is.null(sam.file) & !used.fastq){
    try(system(paste0("awk 'NR > ",sam.indices[4]," {print $",sam.indices[1],"}' ",getwd(),"/",sam.file," > ",getwd(),"/",seq.file)))
    try(system(paste0("awk 'NR > ",sam.indices[4]," {print $",sam.indices[2],"}' ",getwd(),"/",sam.file," > ",getwd(),"/",melo.file)))
    try(system(paste0("awk 'NR > ",sam.indices[4]," {print $",sam.indices[3],"}' ",getwd(),"/",sam.file," > ",getwd(),"/",meca.file)))
  }

  # time stamp for beginning of alignment job
  show(paste0('Align Start:  ',Sys.time()))
  
  # load pre-processed data files
  raw=read.csv(file = seq.file,header = F,sep = ',') # load sequences from pre-processed file
  raw.2=as.matrix(read.csv(melo.file,header = F,sep = ";")[,]) # load CpG indices from pre-processed file
  raw.3=read.csv(meca.file,header = F,sep = ";") # load methyl and hydroxy-methyl scores from pre-processed file

  # convert input reference sequences into vector of bases
  FWD=str_extract_all(parent,boundary("character"))[[1]]
  REV=str_extract_all(target,boundary("character"))[[1]]

  # set read length filters
  if(is.null(read.length)){
    r98=(nchar(raw$V1)[order(nchar(raw$V1))])[round(0.98*length(raw$V1))]
    read.length=c(round(length(REV)/2),r98)
  }

  # load CpG indices, hydroxy-methyl scores, and methyl scores into alignable arrays
  Chm=str_split(raw.2[,1],',')
  Cm=str_split(raw.2[,2],',')
  C.key=str_split(raw.3[,1],',')

  # loop to perform parallelized alignments with scoring on a per-read basis
  registerDoParallel(detectCores())
  if(exact.search==T){
    DATA <- foreach(seq = 1:nrow(raw)) %dopar% {
      find.motif(read=raw[seq,1],Cm = cumsum(as.numeric(Cm[[seq]][-1])+1),Chm = cumsum(as.numeric(Chm[[seq]][-1])+1),round(as.numeric(C.key[[seq]][-1])/256,2),read.length=read.length,FWD=FWD,REV=REV)
    }
  }else{
    if(exact.search==F){
      DATA <- foreach(seq = 1:nrow(raw)) %dopar% {
        map.fragments(read=raw[seq,1],Cm = cumsum(as.numeric(Cm[[seq]][-1])+1),Chm = cumsum(as.numeric(Chm[[seq]][-1])+1),round(as.numeric(C.key[[seq]][-1])/256,2),read.length=read.length,FWD=FWD,REV=REV)
      }
    }
    if(is.integer(exact.search) | is.numeric(exact.search)){
      # prepare kmer sets
      fwd.kset=rep(NA,times=nchar(parent)-exact.search+1)
      rev.kset=rep(NA,times=nchar(target)-exact.search+1)
      for(i in 1:(nchar(parent)-exact.search+1)){
        fwd.kset[i]=paste0(FWD[i:(i+exact.search-1)],collapse = '')
        rev.kset[i]=paste0(REV[i:(i+exact.search-1)],collapse = '')
      }
      share.set=fwd.kset[fwd.kset %in% rev.kset]
      junc.parent=matrix(c(rep(FWD,times=3),rep(REV,times=3),REV,FWD),ncol = 2*length(FWD),byrow = T)
      junction.set=matrix(NA,ncol=exact.search+1,nrow=4)
      for(i in 1:ncol(junction.set)){
        junction.set[,i]=apply(junc.parent[,1:5 + length(FWD) - exact.search - 1 + i],1,paste0,collapse='')
      }
      fwd.kID=which(!(fwd.kset %in% share.set) & !(fwd.kset %in% junction.set) & !(fwd.kset %in% fwd.kset[duplicated(fwd.kset)]))
      rev.kID=which(!(rev.kset %in% share.set) & !(rev.kset %in% junction.set) & !(rev.kset %in% rev.kset[duplicated(rev.kset)]))
      #
      fwd.k=fwd.kset[fwd.kID]
      rev.k=rev.kset[rev.kID]
      if(length(fwd.k)<1 | length(rev.k)<1){
        stop("ERROR! MeLoELF was unable to find at least one unique k-mer each for the FWD/parent and REV/target reference sequences. Try a different k-mer size, or shift to motif search or whole-fragment alignment.")
      }
      keyF <- setNames(fwd.kID,fwd.k)
      keyR <- setNames(rev.kID,rev.k)
      
      # perform kmer-based alignment
      DATA <- foreach(seq = 1:nrow(raw)) %dopar% {
        map.kmers(read=raw[seq,1],Cm = cumsum(as.numeric(Cm[[seq]][-1])+1),Chm = cumsum(as.numeric(Chm[[seq]][-1])+1),round(as.numeric(C.key[[seq]][-1])/256,2),read.length=read.length,FWD=FWD,REV=REV,fwd.k,rev.k,keyF,keyR)
      }
    }
  }

  # get lengths of reads
  lengths.of.reads <- nchar(raw$V1)

  # get fragment numbers
  if(exact.search==T | is.integer(exact.search) | is.numeric(exact.search)){
    READS = rep(NA,times=length(raw$V1))
    READS[lengths.of.reads>=read.length[1] & lengths.of.reads<=read.length[2]] <- as.numeric(do.call(args = lapply(DATA[lengths.of.reads>=read.length[1] & lengths.of.reads<=read.length[2]],'[[','Nfrag'),what = 'c'))
  }else{
    READS <- round(nchar(raw$V1)/length(FWD)+0.4)+padding
  }

  # add final useful data to end of stable variable
  DATA[['N']]=sum(as.numeric(READS[lengths.of.reads>=read.length[1] & lengths.of.reads<=read.length[2]]))+sum(READS==0,na.rm = T)
  DATA[['FWD']]=FWD
  DATA[['REV']]=REV
  DATA[['RLs']]=as.numeric(lengths.of.reads)
  DATA[['RSs']]=raw$V1

  # export generated alignment data to RData file
  save(DATA,file = align.file)

  # time-stamp for end of alignment job
  show(paste0('Align End:  ',Sys.time()))

}

if(process){

  #######################
  ## Fragment Data Processing
  #######################

  # time stamp for beginning of processing job
  show(paste0('Processing Start:  ',Sys.time()))

  # loading RData file containing alignment data
  if (!exists('DATA')) {
    load(align.file)
  }else{
    rm(list=c(setdiff(ls(),c(PAR.list,'DATA'))))
  }

  # set read length filters
  if(is.null(read.length)){
    r98=(DATA[['RLs']][order(as.numeric(DATA[['RLs']]))])[round(0.98*length(DATA[['RLs']]))]
    read.length=c(round(length(DATA[['REV']])/2),r98)
  }

  # generate empty matrices to consolidate data
  FWD.align=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']]));colnames(FWD.align)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.align=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']]));colnames(REV.align)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  FWD.Chm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.Chm)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.Chm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.Chm)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  FWD.Cm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.Cm)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.Cm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.Cm)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  FWD.index=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.index)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.index=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.index)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  Q.reads=matrix(NA,nrow=DATA[['N']],ncol = 4); colnames(Q.reads)<-c('comp','match','id','polyN')
  mapped.frac=rep(NA,times=length(DATA)-5)
  mapped.frac2=rep(NA,times=length(DATA)-5)
  if(pre.ligated){
    remapped.reads=as.data.frame(matrix(NA,nrow=length(DATA[['RSs']]),ncol=4)); colnames(remapped.reads) <- c('Length','Strand','Coverage','Acc')
  }
  COUNTER=1

  # loop for pulling alignment data out of large DATA containers and consolidating it
  for (i in 1:(length(DATA)-5)){

    if(DATA[i]=='blank'){
      next
    }
    if(!is.null(DATA[[i]][['Nfrag']])){
      if(DATA[[i]][['Nfrag']]==0){
        FWD.align[COUNTER,]=rep('x',times=length(DATA[['FWD']]))
        REV.align[COUNTER,]=rep('x',times=length(DATA[['REV']]))
        FWD.index[COUNTER,]=rep(0,times=length(DATA[['FWD']]))
        REV.index[COUNTER,]=rep(0,times=length(DATA[['REV']]))
        Q.reads[COUNTER,4]=i
        mapped.frac[i]=0
        mapped.frac2[i]=0
        COUNTER=COUNTER+1
      }
      if(DATA[[i]][['Nfrag']]==1){
        FWD.align[COUNTER,]=DATA[[i]][['FWD.align']]
        REV.align[COUNTER,]=DATA[[i]][['REV.align']]
        FWD.Chm[COUNTER,]=DATA[[i]][['FWD.Chm']]
        REV.Chm[COUNTER,]=DATA[[i]][['REV.Chm']]
        FWD.Cm[COUNTER,]=DATA[[i]][['FWD.Cm']]
        REV.Cm[COUNTER,]=DATA[[i]][['REV.Cm']]
        Q.reads[COUNTER,1:3]=DATA[[i]][['Q']]
        Q.reads[COUNTER,4]=i
        FWD.index[COUNTER,]=DATA[[i]][['FWD.I']]
        REV.index[COUNTER,]=DATA[[i]][['REV.I']]
        mapped.frac[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']],DATA[[i]][['REV.I']]))))/DATA[['RLs']][i]
        mapped.frac2[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),],DATA[[i]][['REV.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),]))))/DATA[['RLs']][i]
        if(pre.ligated){
          remapped.reads$Acc[i]=1
          if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
            remapped.reads$Strand[i]='FWD'
          }
          if(sum(DATA[[i]]$Q[,3]=='REV',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
            remapped.reads$Strand[i]='REV'
          }
          if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T) >0 & sum(DATA[[i]]$Q[,3]=='REV',na.rm = T) > 0){
            remapped.reads$Strand[i]='Het'
          }
        }
        COUNTER=COUNTER+1
      }
      if(DATA[[i]][['Nfrag']]>1){
        FWD.align[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.align']]
        REV.align[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.align']]
        FWD.Chm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.Chm']]
        REV.Chm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.Chm']]
        FWD.Cm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.Cm']]
        REV.Cm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.Cm']]
        Q.reads[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),1:3]=DATA[[i]][['Q']]
        Q.reads[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),4]=i
        FWD.index[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.I']]
        REV.index[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.I']]
        mapped.frac[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']],DATA[[i]][['REV.I']]))))/DATA[['RLs']][i]
        mapped.frac2[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),],DATA[[i]][['REV.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),]))))/DATA[['RLs']][i]
        if(pre.ligated){
          remapped.reads$Acc[i]=sum(c(DATA[[i]][['FWD.align']],DATA[[i]][['REV.align']]) == c(matrix(c(DATA[['FWD']],DATA[['REV']]),nrow = nrow(DATA[[i]][['REV.align']]),ncol = length(c(DATA[['FWD']],DATA[['REV']])),byrow = T)),na.rm = T)/sum(na.omit(c(DATA[[i]][['FWD.align']],DATA[[i]][['REV.align']])) %in% c('A','C','T','G'))
          if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
            remapped.reads$Strand[i]='FWD'
          }
          if(sum(DATA[[i]]$Q[,3]=='REV',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
            remapped.reads$Strand[i]='REV'
          }
          if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T) >0 & sum(DATA[[i]]$Q[,3]=='REV',na.rm = T) > 0){
            remapped.reads$Strand[i]='Het'
          }
        }
        COUNTER=COUNTER+nrow(DATA[[i]]$Q)
      }
    }
    if(is.null(DATA[[i]][['Nfrag']])){
      FWD.align[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.align']]
      REV.align[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.align']]
      FWD.Chm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.Chm']]
      REV.Chm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.Chm']]
      FWD.Cm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.Cm']]
      REV.Cm[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.Cm']]
      Q.reads[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),1:3]=DATA[[i]][['Q']]
      Q.reads[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),4]=i
      FWD.index[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['FWD.I']]
      REV.index[COUNTER:(COUNTER+length(DATA[[i]][['Q']][,3])-1),]=DATA[[i]][['REV.I']]
      mapped.frac[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']],DATA[[i]][['REV.I']]))))/DATA[['RLs']][i]
      mapped.frac2[i]=sum(!is.na(unique(c(DATA[[i]][['FWD.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),],DATA[[i]][['REV.I']][which(DATA[[i]][['Q']][,1]>=completeness & DATA[[i]][['Q']][,2]>=matching),]))))/DATA[['RLs']][i]
      if(pre.ligated){
        remapped.reads$Acc[i]=sum(c(DATA[[i]][['FWD.align']],DATA[[i]][['REV.align']]) == c(matrix(c(DATA[['FWD']],DATA[['REV']]),nrow = nrow(DATA[[i]][['REV.align']]),ncol = length(c(DATA[['FWD']],DATA[['REV']])),byrow = T)),na.rm = T)/sum(na.omit(c(DATA[[i]][['FWD.align']],DATA[[i]][['REV.align']])) %in% c('A','C','T','G'))
        if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
          remapped.reads$Strand[i]='FWD'
        }
        if(sum(DATA[[i]]$Q[,3]=='REV',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
          remapped.reads$Strand[i]='REV'
        }
        if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T) >0 & sum(DATA[[i]]$Q[,3]=='REV',na.rm = T) > 0){
          remapped.reads$Strand[i]='Het'
        }
      }
      COUNTER=COUNTER+nrow(DATA[[i]]$Q)
    }

  }
  if(pre.ligated){
    remapped.reads$Length[DATA[['RLs']]>=min(read.length) & DATA[['RLs']]<=max(read.length)]=DATA[['RLs']][DATA[['RLs']]>=min(read.length) & DATA[['RLs']]<=max(read.length)]
    remapped.reads$Coverage=mapped.frac
  }

  # prune empty fragment rows
  junk=which(is.na(Q.reads[,3]))
  #
  FWD.index=FWD.index[-junk,]
  REV.index=REV.index[-junk,]
  Q.reads=Q.reads[-junk,]
  FWD.align=FWD.align[-junk,]
  REV.align=REV.align[-junk,]
  FWD.Chm=FWD.Chm[-junk,]
  REV.Chm=REV.Chm[-junk,]
  FWD.Cm=FWD.Cm[-junk,]
  REV.Cm=REV.Cm[-junk,]

  # calculate some QC
  FragPos.fwd=pos.score(FWD.index,DATA[['RLs']][as.numeric(Q.reads[,4])],length(DATA[['FWD']]));colnames(FragPos.fwd)<-c('rloc','edge','5p','3p','mid')
  FragPos.rev=pos.score(REV.index,DATA[['RLs']][as.numeric(Q.reads[,4])],length(DATA[['REV']]));colnames(FragPos.rev)<-c('rloc','edge','5p','3p','mid')
  FragPos.fwdQ=FragPos.fwd[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='FWD'),]
  FragPos.revQ=FragPos.rev[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='REV'),]
  #
  pMAP=sum(mapped.frac*DATA[['RLs']],na.rm = T)/sum(DATA[['RLs']][!is.na(mapped.frac)])
  pMAP2=sum(mapped.frac2*DATA[['RLs']],na.rm = T)/sum(DATA[['RLs']][!is.na(mapped.frac2)])
  pMAP0=sum(DATA[['RLs']][DATA[['RLs']]>=read.length[1] & DATA[['RLs']]<=read.length[2]])/sum(DATA[['RLs']])
  #
  FWD.mat=matrix(DATA[['FWD']],nrow = nrow(FWD.align),ncol = ncol(FWD.align),byrow = T)
  REV.mat=matrix(DATA[['REV']],nrow = nrow(REV.align),ncol = ncol(REV.align),byrow = T)
  FWD.good=colMeans(FWD.align[Q.reads[,3]=='FWD',]==FWD.mat[Q.reads[,3]=='FWD',],na.rm = T)
  FWD.miss=colMeans(FWD.align[Q.reads[,3]=='FWD',]=='x',na.rm = T)
  FWD.err=colMeans(FWD.align[Q.reads[,3]=='FWD',]!=FWD.mat[Q.reads[,3]=='FWD',] & FWD.align[Q.reads[,3]=='FWD',] %in% c('A','C','G','T','N'),na.rm = T)
  REV.good=colMeans(REV.align[Q.reads[,3]=='REV',]==REV.mat[Q.reads[,3]=='REV',],na.rm = T)
  REV.miss=colMeans(REV.align[Q.reads[,3]=='REV',]=='x',na.rm = T)
  REV.err=colMeans(REV.align[Q.reads[,3]=='REV',]!=REV.mat[Q.reads[,3]=='REV',] & REV.align[Q.reads[,3]=='REV',] %in% c('A','C','G','T','N'),na.rm = T)
  rm(FWD.mat,REV.mat)
  #
  FRprop=matrix(table(data.frame(Q.reads[,3],as.numeric(Q.reads[,4]))),nrow = 2)
  FRskew=(FRprop[1,]-FRprop[2,])/colSums(FRprop)
  FRpropQ=matrix(table(data.frame(Q.reads[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),3],as.numeric(Q.reads[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),4]))),nrow = 2)
  FRskewQ=(FRpropQ[1,]-FRpropQ[2,])/colSums(FRpropQ)
  #
  if(methyl.type=='B'){
    Mscore=rowMeans(cbind(FWD.Chm[,FWD.sites]+FWD.Cm[,FWD.sites],REV.Chm[,REV.sites]+REV.Cm[,REV.sites]),na.rm = T)
    MscoreQ=rowMeans(cbind(FWD.Chm[,FWD.sites]+FWD.Cm[,FWD.sites],REV.Chm[,REV.sites]+REV.Cm[,REV.sites]),na.rm = T)[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching)]
  }
  if(methyl.type=='M'){
    Mscore=rowMeans(cbind(FWD.Cm[,FWD.sites],REV.Cm[,REV.sites]),na.rm = T)
    MscoreQ=rowMeans(cbind(FWD.Cm[,FWD.sites],REV.Cm[,REV.sites]),na.rm = T)[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching)]
  }
  if(methyl.type=='H'){
    Mscore=rowMeans(cbind(FWD.Chm[,FWD.sites],REV.Chm[,REV.sites]),na.rm = T)
    MscoreQ=rowMeans(cbind(FWD.Chm[,FWD.sites],REV.Chm[,REV.sites]),na.rm = T)[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching)]
  }
  Mframe=data.frame('S'=Mscore,'N'=as.numeric(Q.reads[,4]))
  MskewZ=aggregate(S ~ N,Mframe,mean)
  Mskew=MskewZ$S
  MframeQ=data.frame('S'=MscoreQ,'N'=as.numeric(Q.reads[which(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),4]))
  MskewZQ=aggregate(S ~ N,MframeQ,mean)
  MskewQ=MskewZQ$S

  if(pre.ligated){

    # compile data for REV polymers with sufficient mapping quality
    polymer.ids=which(remapped.reads$Strand=='REV' & remapped.reads$Coverage>=completeness & remapped.reads$Acc>=matching)
    polymer.Chm=matrix(NA,nrow=length(polymer.ids),ncol = round((max(remapped.reads$Length[polymer.ids],na.rm = T)/length(DATA[['REV']])+1)*length(REV.sites)))
    polymer.Cm=matrix(NA,nrow=length(polymer.ids),ncol = round((max(remapped.reads$Length[polymer.ids],na.rm = T)/length(DATA[['REV']])+1)*length(REV.sites)))
    for (i in 1:length(polymer.ids)){
      temp.Chm=rep(NA,times=remapped.reads$Length[polymer.ids[i]])
      temp.Cm=rep(NA,times=remapped.reads$Length[polymer.ids[i]])
      row.ids=which(as.numeric(Q.reads[,4])==polymer.ids[i])
      temp.Chm[na.omit(c(REV.index[row.ids,REV.sites]))]=c(REV.Chm[row.ids,REV.sites])[!is.na(c(REV.index[row.ids,REV.sites]))]
      temp.Cm[na.omit(c(REV.index[row.ids,REV.sites]))]=c(REV.Cm[row.ids,REV.sites])[!is.na(c(REV.index[row.ids,REV.sites]))]
      temp.Chm=as.numeric(na.omit(temp.Chm))
      temp.Cm=as.numeric(na.omit(temp.Cm))
      polymer.Chm[i,1:length(temp.Chm)]=temp.Chm
      polymer.Cm[i,1:length(temp.Cm)]=temp.Cm
    }
    polymer.Both=polymer.Chm+polymer.Cm

    # compile data for FWD polymers with sufficient mapping quality
    polymer.ids.fwd=which(remapped.reads$Strand=='FWD' & remapped.reads$Coverage>=completeness & remapped.reads$Acc>=matching)
    polymer.Chm.fwd=matrix(NA,nrow=length(polymer.ids.fwd),ncol = round((max(remapped.reads$Length[polymer.ids.fwd],na.rm = T)/length(DATA[['FWD']])+1)*length(FWD.sites)))
    polymer.Cm.fwd=matrix(NA,nrow=length(polymer.ids.fwd),ncol = round((max(remapped.reads$Length[polymer.ids.fwd],na.rm = T)/length(DATA[['FWD']])+1)*length(FWD.sites)))
    for (i in 1:length(polymer.ids.fwd)){
      temp.Chm=rep(NA,times=remapped.reads$Length[polymer.ids.fwd[i]])
      temp.Cm=rep(NA,times=remapped.reads$Length[polymer.ids.fwd[i]])
      row.ids.fwd=which(as.numeric(Q.reads[,4])==polymer.ids.fwd[i])
      temp.Chm[na.omit(c(FWD.index[row.ids.fwd,FWD.sites]))]=c(FWD.Chm[row.ids.fwd,FWD.sites])[!is.na(c(FWD.index[row.ids.fwd,FWD.sites]))]
      temp.Cm[na.omit(c(FWD.index[row.ids.fwd,FWD.sites]))]=c(FWD.Cm[row.ids.fwd,FWD.sites])[!is.na(c(FWD.index[row.ids.fwd,FWD.sites]))]
      temp.Chm=as.numeric(na.omit(temp.Chm))
      temp.Cm=as.numeric(na.omit(temp.Cm))
      polymer.Chm.fwd[i,1:length(temp.Chm)]=temp.Chm
      polymer.Cm.fwd[i,1:length(temp.Cm)]=temp.Cm
    }
    polymer.Both.fwd=polymer.Chm.fwd+polymer.Cm.fwd

    # clean up data sets
    if(methyl.type=='B'){
      polymer.actual.fwd=polymer.Both.fwd
      polymer.actual.rev=polymer.Both
    }
    if(methyl.type=='M'){
      polymer.actual.fwd=polymer.Cm.fwd
      polymer.actual.rev=polymer.Cm
    }
    if(methyl.type=='H'){
      polymer.actual.fwd=polymer.Chm.fwd
      polymer.actual.rev=polymer.Chm
    }
    #
    save(polymer.actual.rev,polymer.actual.fwd,file = processed.file)

  }

  # pull methylation data for only reads with sufficient mapping quality
  if(!pre.ligated){
    FWD.Chm.pruned=FWD.Chm[which((Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='FWD')),]
    FWD.Cm.pruned=FWD.Cm[which((Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='FWD')),]
    REV.Chm.pruned=REV.Chm[which((Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='REV')),]
    REV.Cm.pruned=REV.Cm[which((Q.reads[,1]>=completeness & Q.reads[,2]>=matching & Q.reads[,3]=='REV')),]
  }

  # calculate total methylation data
  if(!pre.ligated){
    FWD.both.pruned=FWD.Chm.pruned+FWD.Cm.pruned
    try(FWD.both.pruned[is.na(FWD.Cm.pruned) & !is.na(FWD.Chm.pruned)]<-FWD.Chm.pruned[is.na(FWD.Cm.pruned) & !is.na(FWD.Chm.pruned)])
    try(FWD.both.pruned[!is.na(FWD.Cm.pruned) & is.na(FWD.Chm.pruned)]<-FWD.Cm.pruned[!is.na(FWD.Cm.pruned) & is.na(FWD.Chm.pruned)])
    REV.both.pruned=REV.Chm.pruned+REV.Cm.pruned
    try(REV.both.pruned[is.na(REV.Cm.pruned) & !is.na(REV.Chm.pruned)]<-REV.Chm.pruned[is.na(REV.Cm.pruned) & !is.na(REV.Chm.pruned)])
    try(REV.both.pruned[!is.na(REV.Cm.pruned) & is.na(REV.Chm.pruned)]<-REV.Cm.pruned[!is.na(REV.Cm.pruned) & is.na(REV.Chm.pruned)])
  }

  # clean up data sets
  if(!pre.ligated){
    if(methyl.type=='B'){
      data.actual.rev=REV.both.pruned[rowSums(!is.na(REV.both.pruned))>0,rev(REV.sites)]
      data.actual.fwd=FWD.both.pruned[rowSums(!is.na(FWD.both.pruned))>0,FWD.sites]
    }
    if(methyl.type=='M'){
      data.actual.rev=REV.Cm.pruned[rowSums(!is.na(REV.Cm.pruned))>0,rev(REV.sites)]
      data.actual.fwd=FWD.Cm.pruned[rowSums(!is.na(FWD.Cm.pruned))>0,FWD.sites]
    }
    if(methyl.type=='H'){
      data.actual.rev=REV.Chm.pruned[rowSums(!is.na(REV.Chm.pruned))>0,rev(REV.sites)]
      data.actual.fwd=FWD.Chm.pruned[rowSums(!is.na(FWD.Chm.pruned))>0,FWD.sites]
    }
    #
    save(data.actual.fwd,data.actual.rev,file = processed.file)
  }

  #######################
  ## Processivity Analysis and Graphing
  #######################

  if(pre.ligated){
    # Get threshold values
    if(thresh.meth=='AUC'){
      thresh <- find_thresh_for_auc(polymer.actual.fwd, polymer.actual.rev)
    }
    if(thresh.meth=='BM'){
      thresh=BM.thresh(polymer.actual.fwd,polymer.actual.rev)
    }
    if(is.numeric(thresh.meth)){
      thresh=thresh.meth
    }

    fwd.binary <- binarize_matrix(polymer.actual.fwd, thresh)
    rev.binary <- binarize_matrix(polymer.actual.rev, thresh)
  }
  if(!pre.ligated){
    # Get threshold values
    if(thresh.meth=='AUC'){
      thresh <- find_thresh_for_auc(data.actual.fwd, data.actual.rev)
    }
    if(thresh.meth=='BM'){
      if(var.thresh<2){
        thresh=BM.thresh(data.actual.fwd,data.actual.rev)
      }
      if(var.thresh>0){
        threshIrev=rep(NA,times=ncol(data.actual.rev))
        if(var.thresh==2){
          threshIfwd=threshIrev
        }
        for(i in 1:length(threshIrev)){
          threshIrev[i]=BM.thresh(data.actual.fwd[,i],data.actual.rev[,i],set = 'r')
          if(var.thresh==2){
            threshIfwd[i]=BM.thresh(data.actual.fwd[,i],data.actual.rev[,i],set = 'f')
          }
        }
      }
    }
    if(is.numeric(thresh.meth)){
      thresh=thresh.meth
    }

    if(var.thresh==0 | thresh.meth=='AUC' | is.numeric(thresh.meth)){
      fwd.binary <- binarize_matrix(data.actual.fwd, thresh)
      rev.binary <- binarize_matrix(data.actual.rev, thresh)
    }else{
      if(var.thresh==2){
        fwd.binary <- binarize_matrix(data.actual.fwd, threshIfwd)
      }else{
        fwd.binary <- binarize_matrix(data.actual.fwd, thresh)
      }
      rev.binary <- binarize_matrix(data.actual.rev, threshIrev)
    }
  }
  ctrl.binary=rev.binary; ctrl.binary[!is.na(ctrl.binary)]=0
  ctrl.binary[sample(which(!is.na(rev.binary)),sum(rev.binary==1,na.rm = T))]=1

  fwd.binary.surv <- get_survival_data(fwd.binary)
  rev.binary.surv <- get_survival_data(rev.binary[,ncol(rev.binary):1])
  ctrl.binary.surv <- get_survival_data(ctrl.binary[,ncol(ctrl.binary):1])

  fwd.surv.plot <- empirical_survival(fwd.binary.surv$duration)
  rev.surv.plot <- empirical_survival(rev.binary.surv$duration)
  ctrl.surv.plot <- empirical_survival(ctrl.binary.surv$duration)

  fwd.frac.full <- fraction_full(fwd.binary.surv)
  rev.frac.full <- fraction_full(rev.binary.surv)
  ctrl.frac.full <- fraction_full(ctrl.binary.surv)
  frac.full <- c(fwd.frac.full, rev.frac.full, ctrl.frac.full)

  fwd.auc <- compute_auc(fwd.surv.plot)
  rev.auc <- compute_auc(rev.surv.plot)
  ctrl.auc <- compute_auc(ctrl.surv.plot)
  auc <- c(fwd.auc, rev.auc, ctrl.auc)

  fwd.median.dur <- median_duration(fwd.surv.plot)
  rev.median.dur <- median_duration(rev.surv.plot)
  ctrl.median.dur <- median_duration(ctrl.surv.plot)
  median <- c(fwd.median.dur, rev.median.dur, ctrl.median.dur)

  matrices <- list(FWD = fwd.binary.surv, REV = rev.binary.surv, CTRL = ctrl.binary.surv)
  summary<- data.frame(strand = names(matrices), frac.full, auc, median)

  save(matrices, file="methylation_matrices.RDS")
  write.csv(summary, "survival_summary.csv")

  # calculate Modkit-like overall methylation for reads
  check.melo=fread(file = melo.file,header = F,sep = ';',fill = T)[,-3]
  check.meca=fread(file = meca.file,header = F,sep = ';')[,1]
  check.seq=fread(file = seq.file,header = F,sep = ';')[,1]
  #
  registerDoParallel(detectCores())
  if(var.thresh<2 | pre.ligated){
    modk.dat <- foreach(i = 1:length(check.seq),.combine = 'rbind') %dopar% {
      modk.sum(melo = check.melo[i,],meca = check.meca[i],seqq = check.seq[i],methyl.type = methyl.type,thresh = thresh,fill.Cs = modk.fillC)
    }
  }
  if(var.thresh==2 & !(pre.ligated)){
    modk.dat <- foreach(i = 1:length(check.seq),.combine = 'rbind') %dopar% {
      modk.sum(melo = check.melo[i,],meca = check.meca[i],seqq = check.seq[i],methyl.type = methyl.type,thresh = mean(threshIrev),fill.Cs = modk.fillC)
    }
  }
  rm(check.meca,check.melo,check.seq)
  #
  used.set=which(DATA[['RLs']]>=min(read.length) & DATA[['RLs']]<=max(read.length))
  all.modk=sum(modk.dat$pC*modk.dat$Nc,na.rm = T)/sum(modk.dat$Nc)
  used.modk=sum(modk.dat$pC[used.set]*modk.dat$Nc[used.set],na.rm = T)/sum(modk.dat$Nc[used.set])
  filt.modk=sum(modk.dat$pC[-used.set]*modk.dat$Nc[-used.set],na.rm = T)/sum(modk.dat$Nc[-used.set])
  rm(modk.dat)

  # quality control analyses
  if(!pre.ligated){
    read.filt=nrow(data.actual.rev)/sum(Q.reads[,3]=='REV',na.rm = T)
    read.filt.fwd=nrow(data.actual.fwd)/sum(Q.reads[,3]=='FWD',na.rm = T)
    fNA.fwd=colSums(is.na(fwd.binary))/nrow(fwd.binary)
    fNA.rev=colSums(is.na(rev.binary))/nrow(rev.binary)
    fCpGs.fwd=colSums(fwd.binary,na.rm = T)/colSums(!is.na(fwd.binary))
    fCpGs.rev=colSums(rev.binary,na.rm = T)/colSums(!is.na(rev.binary))
  }
  fSs.fwd=sum(rowSums(fwd.binary,na.rm = T)>0)/nrow(fwd.binary)
  fSs.rev=sum(rowSums(rev.binary,na.rm = T)>0)/nrow(rev.binary)
  fMs.fwd=sum(fwd.binary,na.rm = T)/sum(!is.na(fwd.binary))
  fMs.rev=sum(rev.binary,na.rm = T)/sum(!is.na(rev.binary))
  if(!pre.ligated){
    f53m.rev=rep(NA,times=nrow(rev.binary))
    f35m.rev=rep(NA,times=nrow(rev.binary))
    for(i in 1:nrow(rev.binary)){
      if(sum(rev.binary[i,]==1,na.rm = T)>0){
        f53m.rev[i]=max(which(rev.binary[i,]==1))
        f35m.rev[i]=min(which(rev.binary[i,]==1))
      }
    }
    f53m.rev=table(f53m.rev)/sum(!is.na(f53m.rev))
    temp=rep(0,times=length(REV.sites));temp[as.numeric(names(f53m.rev))]=as.numeric(f53m.rev)
    f53m.rev=temp
    f35m.rev=table(f35m.rev)/sum(!is.na(f35m.rev))
    temp=rep(0,times=length(REV.sites));temp[as.numeric(names(f35m.rev))]=as.numeric(f35m.rev)
    f35m.rev=temp
  }
  #
  if(!pre.ligated){
    png('MethylScoring.png', height = 2650, width = 2800, res=300)
    par(mfrow=c(2,1),mar=c(4,5,3,1))
    #
    FWD.Cm.pdfs=list(NULL)
    REV.Cm.pdfs=list(NULL)
    for (i in 1:length(FWD.sites)){
      FWD.Cm.pdfs[[i]]=density(data.actual.fwd[,i],na.rm = T,from = 0,to = 1)
      REV.Cm.pdfs[[i]]=density(data.actual.rev[,i],na.rm = T,from = 0,to = 1)
    }
    plot(NULL,NULL,ylim=c(0,1),xlim=c(0,length(FWD.sites)),main = paste0(plot_title,'Methyl Score Distributions'),cex.axis = 1.6,ylab = 'Methyl Score',cex.lab=2,cex.main=2,xaxt='n',xlab='')
    if(var.thresh<2 | thresh.meth=='AUC'){
      abline(h=thresh,lty='dashed',col='grey',lwd=2)
    }
    for(i in 1:length(FWD.sites)){
      lines(i-0.725-FWD.Cm.pdfs[[i]][['y']]/max(FWD.Cm.pdfs[[i]][['y']])*0.2,FWD.Cm.pdfs[[i]][['x']],col='blue')
      lines(i-0.725+FWD.Cm.pdfs[[i]][['y']]/max(FWD.Cm.pdfs[[i]][['y']])*0.2,FWD.Cm.pdfs[[i]][['x']],col='blue')
      lines(i-0.275-REV.Cm.pdfs[[i]][['y']]/max(REV.Cm.pdfs[[i]][['y']])*0.2,REV.Cm.pdfs[[i]][['x']],col='red')
      lines(i-0.275+REV.Cm.pdfs[[i]][['y']]/max(REV.Cm.pdfs[[i]][['y']])*0.2,REV.Cm.pdfs[[i]][['x']],col='red')
      if(var.thresh>0 & thresh.meth=='BM'){
        lines(i-0.275+c(-1/5,1/5),rep(threshIrev[i],times=2),lwd=4,col='grey')
        if(var.thresh==2){
          lines(i-0.725+c(-1/5,1/5),rep(threshIfwd[i],times=2),lwd=4,col='grey')
        }
      }
    }
    points(1:length(FWD.sites)-0.725,apply(data.actual.fwd,2,median,na.rm=T),col='purple',pch='-',cex=5)
    points(1:length(REV.sites)-0.275,apply(data.actual.rev,2,median,na.rm=T),col='purple',pch='-',cex=5)
    axis(side = 1,at = c(1:length(FWD.sites))-0.5,labels = paste0(DATA[['FWD']][FWD.sites],'p',DATA[['FWD']][FWD.sites+1],'-',FWD.sites+0.5),cex.axis=1.2)
    axis(side = 1,at = c(1:length(FWD.sites))-0.725,labels = paste0(round(colSums(!is.na(data.actual.fwd))/1e3),'k'),tick = F,line = 1.2,cex.axis=1.0,col.axis = 'blue')
    axis(side = 1,at = c(1:length(FWD.sites))-0.275,labels = paste0(round(colSums(!is.na(data.actual.rev))/1e3),'k'),tick = F,line = 1.2,cex.axis=1.0,col.axis = 'red')
    legend('topright',legend = plot.nom,col = c('blue','red'),fill = c('blue','red'),cex=1.1,bty = 'n')
    #
    par(mar=c(3,5,3,1))
    plot(NULL,NULL,xlim=c(0,5*length(FWD.sites)+4),ylim=c(0,1.1),xaxt='n',ylab='Fraction of Product',main=paste0(plot_title,'Methyl Location Distributions'),cex.main=2,cex.lab=1.5,cex.axis=1.5,xlab='')
    axis(side = 1,at = c(0:6)*5+2,labels = paste0(DATA[['FWD']][FWD.sites],'p',DATA[['FWD']][FWD.sites+1],'-',FWD.sites+0.5),cex.axis=1.1)
    abline(h=fMs.rev,col='red',lty='dashed',lwd=2)
    abline(h=fSs.rev,col='purple',lty='dashed',lwd=2)
    if(!(exact.search==T)){
      abline(h=read.filt,col='orange',lty='dotted',lwd=2)
    }
    abline(h=pMAP,col='cyan',lty='dotted',lwd=2)
    points(c(0:6)*5+0.5,fCpGs.fwd,type='h',pch=22,lwd=15,col=1)
    points(c(0:6)*5+1.5,fCpGs.rev,type='h',pch=22,lwd=15,col=2)
    points(c(0:6)*5+2.5,f53m.rev,type='h',pch=22,lwd=15,col=3)
    points(c(0:6)*5+3.5,f35m.rev,type='h',pch=22,lwd=15,col=4)
    legend('topright',legend=c(paste0(plot.nom,' Methyls'),"5'->3' Start","3'->5' Start"),col=1:4,fill=1:4,cex=0.8,bty = 'n')
    if(exact.search==T){
      text(x=(5*length(FWD.sites)+4)*1.04,y=c(0.8,0.7),pos = 2,col = c('red','purple'),cex = 1.0,labels = paste0('(',round(100*c(fMs.fwd,fSs.fwd)),'%) ',round(100*c(fMs.rev,fSs.rev)),c('% CpG','% Sub.')))
    }else{
      text(x=(5*length(FWD.sites)+4)*1.04,y=c(0.8,0.7,0.6),pos = 2,col = c('red','purple','orange'),cex = 1.0,labels = paste0('(',round(100*c(fMs.fwd,fSs.fwd,read.filt.fwd)),'%) ',round(100*c(fMs.rev,fSs.rev,read.filt)),c('% CpG','% Sub.','% Qual.')))
    }
    text(x=0,y=1.07,pos=4,labels=paste0('[',round(100*pMAP0),'%]  ',round(100*pMAP),'/',round(100*pMAP2),'% Map'),col='cyan',cex=1.3)
    text(x=2*length(FWD.sites),y=1.07,pos=4,labels=paste0('[',round(100*filt.modk),'/',round(100*all.modk),'%]  ',round(100*used.modk),'% 5(h)mC'),col='grey',cex=1.3)
    #
    dev.off()
  }
  png('FragPos.png', height = 2650, width = 2800, res=300)
  #
  par(fig=c(0,1,0.67,1),mar=c(5,5,3,1))
  f=1
  plot(pdf.make(c(FragPos.fwdQ[,f],FragPos.revQ[,f])),ylim=c(0,max(c(pdf.make(FragPos.fwdQ[,f])$y,pdf.make(FragPos.revQ[,f])$y))),xlim=c(0,1),lwd=3,xlab="Read Location (5' -> 3')",ylab='Density',main = paste0(plot_title,'Fragment Positions: Relative to Read'))
  lines(pdf.make(FragPos.fwdQ[,f]),col='blue')
  lines(pdf.make(FragPos.revQ[,f]),col='red')
  legend('topright',legend = c('Both',plot.nom),fill = c('black','blue','red'),col=c('black','blue','red'),bty = 'n')
  #
  par(fig=c(0,0.5,0.33,0.67),mar=c(5,5,3,1),new=T)
  f=2
  plot(pdf.make(c(FragPos.fwdQ[,f],FragPos.revQ[,f])),ylim=c(0,max(c(pdf.make(FragPos.fwdQ[,f])$y,pdf.make(FragPos.revQ[,f])$y))),xlim=c(0,5),lwd=3,xlab='Fragments Away',ylab='Density',main = 'From Nearest Edge')
  lines(pdf.make(FragPos.fwdQ[,f]),col='blue')
  lines(pdf.make(FragPos.revQ[,f]),col='red')
  #
  par(fig=c(0.5,1,0.33,0.67),mar=c(5,5,3,1),new=T)
  f=3
  plot(pdf.make(c(FragPos.fwdQ[,f],FragPos.revQ[,f])),ylim=c(0,max(c(pdf.make(FragPos.fwdQ[,f])$y,pdf.make(FragPos.revQ[,f])$y))),xlim=c(0,10),lwd=3,xlab='Fragments Away',ylab='Density',main = "From 5' End")
  lines(pdf.make(FragPos.fwdQ[,f]),col='blue')
  lines(pdf.make(FragPos.revQ[,f]),col='red')
  #
  par(fig=c(0,0.5,0,0.33),mar=c(5,5,3,1),new=T)
  f=4
  plot(pdf.make(c(FragPos.fwdQ[,f],FragPos.revQ[,f])),ylim=c(0,max(c(pdf.make(FragPos.fwdQ[,f])$y,pdf.make(FragPos.revQ[,f])$y))),xlim=c(-10,0),lwd=3,xlab='Fragments Away',ylab='Density',main = "From 3' End")
  lines(pdf.make(FragPos.fwdQ[,f]),col='blue')
  lines(pdf.make(FragPos.revQ[,f]),col='red')
  #
  par(fig=c(0.5,1,0,0.33),mar=c(5,5,3,1),new=T)
  f=5
  plot(pdf.make(c(FragPos.fwdQ[,f],FragPos.revQ[,f])),ylim=c(0,max(c(pdf.make(FragPos.fwdQ[,f])$y,pdf.make(FragPos.revQ[,f])$y))),xlim=c(-5,5),lwd=3,xlab='Fragments Away',ylab='Density',main = "From Read Center")
  lines(pdf.make(FragPos.fwdQ[,f]),col='blue')
  lines(pdf.make(FragPos.revQ[,f]),col='red')
  #
  dev.off()
  #
  if(!(exact.search==T)){
    png('FragSeq.png', height = round(2650*0.7), width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(6,3,3,1))
    #
    plot(NULL,NULL,ylim=c(0,1),xlim=c(0,length(DATA[['FWD']])+1),main = paste0(plot_title,'Substrate Sequencing Accuracy'),cex.axis = 1.3,ylab = '',yaxt='n',cex.lab=1.3,cex.main=2,xaxt='n',xlab='')
    abline(v=FWD.sites+1/2,lty='dotted',col='grey',lwd=2)
    lines(FWD.good,col='blue',lty='solid',lwd=5)
    lines(FWD.err,col='blue',lty='solid',lwd=5)
    lines(FWD.miss,col='blue',lty='solid',lwd=5)
    lines(rev(REV.good),col='red',lty='solid',lwd=5)
    lines(rev(REV.err),col='red',lty='solid',lwd=5)
    lines(rev(REV.miss),col='red',lty='solid',lwd=5)
    lines(FWD.good,col='grey',lty='solid',lwd=1.5)
    lines(FWD.err,col='black',lty='solid',lwd=1.5)
    lines(FWD.miss,col='white',lty='solid',lwd=1.5)
    lines(rev(REV.good),col='grey',lty='solid',lwd=1.5)
    lines(rev(REV.err),col='black',lty='solid',lwd=1.5)
    lines(rev(REV.miss),col='white',lty='solid',lwd=1.5)
    legend('topright',legend = c(plot.nom,'Correct','Miscall','Deletion'),col = c('blue','red','grey','black','white'),fill=c('blue','red','grey','black','white'),bty = 'n',cex=1)
    axis(side = 2,line = 0,at = c(0,0.25,0.5,3/4,1),labels = c('0%','25%','50%','75%','100%'),tick = T,cex.axis=1.3,col.axis = 'black')
    axis(side = 1,at = c(0:(length(DATA[['FWD']])+1)),labels = c("5'",DATA[['FWD']],"3'"),cex.axis=0.7,col.axis = 'blue')
    axis(side = 1,line = 1,at = c(0:(length(DATA[['REV']])+1)),labels = rev(c("5'",DATA[['REV']],"3'")),tick = F,cex.axis=0.7,col.axis = 'red')
    #
    dev.off()
  }
  #
  if(T){
    png('FwdRevReadBias.png', height = round(2650*0.8), width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(5,5,3,1))
    #
    plot(NULL,NULL,ylim=c(0,max(c(pdf.make(FRskew,pars = c(-1.1,1.1,0.05))$y,pdf.make(FRskewQ,pars = c(-1.1,1.1,0.05))$y))),xlim=c(-1.1,1.1),main = paste0(plot_title,plot.nom[1],' vs ',plot.nom[2],' Per-Read Bias'),cex.axis = 1.4,ylab = 'Probability Density',cex.lab=1.6,cex.main=2,xlab=paste0('Proportion Bias (',plot.nom[2],' <--> ',plot.nom[1],')'))
    lines(pdf.make(FRskew,pars = c(-1.1,1.1,0.05)),col='black',lty='solid',lwd=4)
    lines(pdf.make(FRskewQ,pars = c(-1.1,1.1,0.05)),col='purple',lty='solid',lwd=4)
    legend('top',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    #
    dev.off()
  }
  #
  if(T){
    png('MethylReadBias.png', height = round(2650*0.8), width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(5,5,3,1))
    #
    plot(NULL,NULL,ylim=c(0,max(c(pdf.make(Mskew,pars = c(0,1,0.025))$y,pdf.make(MskewQ,pars = c(0,1,0.025))$y))),xlim=c(0,1),main = paste0(plot_title,'CpG vs 5(h)mCpG Per-Read Bias'),cex.axis = 1.4,ylab = 'Probability Density',cex.lab=1.6,cex.main=2,xlab=paste0('Fragment-CpG Methyl Score Average'))
    lines(pdf.make(Mskew,pars = c(0,1,0.025)),col='black',lty='solid',lwd=4)
    lines(pdf.make(MskewQ,pars = c(0,1,0.025)),col='purple',lty='solid',lwd=4)
    legend('topright',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    #
    dev.off()
  }
  #
  if(T){
    png('FRbiasLengthCorr.png', height = round(2650*0.8), width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(5,5,3,1))
    #
    x1=DATA[['RLs']][unique(na.omit(as.numeric(Q.reads[,4])))]
    x2=DATA[['RLs']][as.numeric(unique(Q.reads[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),4]))]
    #
    plot(NULL,NULL,ylim=c(-1,1),xlim=c(0,max(x1)),main = paste0(plot_title,plot.nom[1],':',plot.nom[2],' Polymer Bias Correlation to Read Length'),cex.axis = 1.4,ylab = paste0('Proportion Bias (',plot.nom[2],' <--> ',plot.nom[1],')'),cex.lab=1.6,cex.main=1.5,xlab=paste0('Read Length (bp)'))
    contour(kde2d(x1,FRskew),col='black',lwd=1.5,add=T)
    contour(kde2d(x2,FRskewQ),col='purple',lwd=1,add=T)
    lines(smooth.spline(x1,FRskew,spar = 0.7),col='black',lty='solid',lwd=4)
    lines(smooth.spline(x2,FRskewQ,spar = 0.7),col='purple',lty='solid',lwd=4)
    legend('topright',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    #
    dev.off()
  }
  #
  if(T){
    png('5mCbiasLengthCorr.png', height = round(2650*0.8), width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(5,5,3,1))
    #
    x1=DATA[['RLs']][as.numeric(MskewZ$N)]
    x2=DATA[['RLs']][as.numeric(MskewZQ$N)]
    #
    plot(NULL,NULL,ylim=c(0,1),xlim=c(0,max(x1)),main = paste0(plot_title,'CpG:5(h)mCpG Bias Correlation to Read Length'),cex.axis = 1.4,ylab = paste0('Fragment-CpG Methyl Score Average'),cex.lab=1.6,cex.main=1.5,xlab=paste0('Read Length (bp)'))
    contour(kde2d(x1,MskewZ$S),col='black',lwd=1.5,add=T)
    contour(kde2d(x2,MskewZQ$S),col='purple',lwd=1,add=T)
    lines(smooth.spline(x1,MskewZ$S,spar = 0.7),col='black',lty='solid',lwd=4)
    lines(smooth.spline(x2,MskewZQ$S,spar = 0.7),col='purple',lty='solid',lwd=4)
    legend('topright',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    #
    dev.off()
  }
  #
  if(T){
    png('5mCReadPosCorr.png', height = round(2650*1.4), width = 2800, res=300)
    par(mfrow=c(2,1),mar=c(5,5,3,1))
    #
    x1=c(FWD.index[,FWD.sites])
    if(methyl.type=='B'){
      y1=c((FWD.Chm+FWD.Cm)[,FWD.sites])
    }
    if(methyl.type=='M'){
      y1=c(FWD.Cm[,FWD.sites])
    }
    if(methyl.type=='H'){
      y1=c(FWD.Chm[,FWD.sites])
    }
    #
    x2=c(FWD.index[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),FWD.sites])
    if(methyl.type=='B'){
      y2=c((FWD.Chm+FWD.Cm)[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),FWD.sites])
    }
    if(methyl.type=='M'){
      y2=c(FWD.Cm[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),FWD.sites])
    }
    if(methyl.type=='H'){
      y2=c(FWD.Chm[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),FWD.sites])
    }
    d1=na.omit(data.frame('x'=x1,'y'=y1))
    d2=na.omit(data.frame('x'=x2,'y'=y2))
    #
    plot(NULL,NULL,xlim=c(0,max(x1,na.rm = T)),ylim=c(0,1),main = paste0(plot_title,'Methylation Correlation to Read Position: ',plot.nom[1]),cex.axis = 1.4,ylab = paste0('Methyl Score'),cex.lab=1.6,cex.main=1.5,xlab=paste0("Bases from 5' End"))
    contour(kde2d(d1$x,d1$y),col='black',lwd=1.5,add=T)
    contour(kde2d(d2$x,d2$y),col='purple',lwd=1,add=T)
    lines(smooth.spline(d1$x,d1$y,spar = 0.7),col='black',lty='solid',lwd=4)
    lines(smooth.spline(d2$x,d2$y,spar = 0.7),col='purple',lty='solid',lwd=4)
    legend('topright',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    #
    #
    x1=c(REV.index[,REV.sites])
    if(methyl.type=='B'){
      y1=c((REV.Chm+REV.Cm)[,REV.sites])
    }
    if(methyl.type=='M'){
      y1=c(REV.Cm[,REV.sites])
    }
    if(methyl.type=='H'){
      y1=c(REV.Chm[,REV.sites])
    }
    #
    x2=c(REV.index[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),REV.sites])
    if(methyl.type=='B'){
      y2=c((REV.Chm+REV.Cm)[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),REV.sites])
    }
    if(methyl.type=='M'){
      y2=c((REV.Cm)[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),REV.sites])
    }
    if(methyl.type=='H'){
      y2=c((REV.Chm)[(Q.reads[,1]>=completeness & Q.reads[,2]>=matching),REV.sites])
    }
    d1=na.omit(data.frame('x'=x1,'y'=y1))
    d2=na.omit(data.frame('x'=x2,'y'=y2))
    #
    plot(NULL,NULL,xlim=c(0,max(x1,na.rm = T)),ylim=c(0,1),main = paste0(plot_title,'Methylation Correlation to Read Position: ',plot.nom[2]),cex.axis = 1.4,ylab = paste0('Methyl Score'),cex.lab=1.6,cex.main=1.5,xlab=paste0("Bases from 5' End"))
    contour(kde2d(d1$x,d1$y),col='black',lwd=1.5,add=T)
    contour(kde2d(d2$x,d2$y),col='purple',lwd=1,add=T)
    lines(smooth.spline(d1$x,d1$y,spar = 0.7),col='black',lty='solid',lwd=4)
    lines(smooth.spline(d2$x,d2$y,spar = 0.7),col='purple',lty='solid',lwd=4)
    legend('topright',legend = c('Pre-Quality Filtering','Post-Quality Filtering'),col = c('black','purple'),fill=c('black','purple'),bty = 'n',cex=1.5)
    dev.off()
  }
  #
  if(pre.ligated){
    png('MethylScoring.png', height = 1325, width = 2800, res=300)
    par(mfrow=c(1,1),mar=c(4,5,3,1))
    #
    FWD.Cm.pdfs=pdf.make(polymer.actual.fwd,pars = c(0,1))
    REV.Cm.pdfs=pdf.make(polymer.actual.rev,pars = c(0,1))
    plot(NULL,NULL,ylim=c(0,max(c(FWD.Cm.pdfs$y,REV.Cm.pdfs$y))),xlim=c(0,1),main = paste0(plot_title,'Methyl Score Distributions'),cex.axis = 1.5,xlab = 'Methyl Score',cex.lab=1.5,cex.main=2,ylab='Probability Density')
    abline(v=thresh,lty='dashed',col='grey',lwd=2)
    lines(FWD.Cm.pdfs,col='blue',lwd=2)
    lines(REV.Cm.pdfs,col='red',lwd=2)
    legend('topright',legend = plot.nom,col = c('blue','red'),fill = c('blue','red'),cex=1.3,bty = 'n')
    #
    dev.off()
  }

  wd <- getwd()

  # plotting survival analysis
  png(paste0(wd, "/MethylationSurvival.png"), height = 1320, width = 2800, res=300)
  par(mfrow=c(1,1),mar=c(5,6,3,8))
  plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,100),
       xlab="Survival score", ylab="Percent of Reads",
       main=paste0(plot_title, "Methylation Survival"), cex.main = 2, cex.axis=2, cex.lab = 2)
  lines(fwd.surv.plot$x, fwd.surv.plot$y, type="s", col="blue", lwd=2)
  lines(rev.surv.plot$x, rev.surv.plot$y, type="s", col="magenta1", lwd=2)
  lines(ctrl.surv.plot$x, ctrl.surv.plot$y, type="s", col="grey", lwd=2)
  legend('topright', legend = c(paste0(plot.nom[1],' (AUC=',round(fwd.auc, digits=2),')'), paste0(plot.nom[2],' (AUC=',round(rev.auc, digits=2),')'), paste0('Sim. Dist. (AUC=',round(ctrl.auc, digits=2),')')),col = c('blue','magenta1','grey'),fill = c('blue','magenta1','grey'),cex=1.2, bty="n")
  dev.off()

  write.csv(rev.binary, "revbinary.csv")

  if(!pre.ligated){
    # Distribution analysis
    fwd_vals <- (colSums(fwd.binary) / nrow(fwd.binary))*100
    rev_vals <- (colSums(rev.binary) / nrow(rev.binary))*100
    ctrl_vals <- (colSums(ctrl.binary) / nrow(ctrl.binary))*100
    n_positions <- length(fwd_vals)

    fwd.binary.counts <- row_sum_counts(na.omit(fwd.binary))
    rev.binary.counts <- row_sum_counts(na.omit(rev.binary))
    ctrl.binary.counts <- row_sum_counts(na.omit(ctrl.binary))

    valsmat <- rbind(fwd_vals, rev_vals,ctrl_vals)
    valsmat_sub <- valsmat[, -1, drop=FALSE]

    new_labels <- 1:(ncol(valsmat) - 1)

    fwd_vals <- fwd.binary.counts[-1, 3]
    rev_vals <- rev.binary.counts[-1, 3]
    ctrl_vals <- ctrl.binary.counts[-1, 3]
    x_labels <- fwd.binary.counts[-1, 1]

    png(paste0(wd, "/MethylDistribution.png"), height = 1320, width = 2800, res=300)
    par(mfrow=c(1,1), mar=c(5,6,3,8), xpd=TRUE)


    suppressWarnings(
    barplot(fwd_vals,ylim = c(0, 100),width = 1,space = c(0, rep(3, times = length(fwd_vals) - 1)),col = 'blue',cex.main = 2,cex.axis = 2,ylab = 'Percent of Reads',xlab = 'Number of 5mCs',cex.lab = 2,main = paste0(plot_title, "Methylation Distribution")))
    barplot(rev_vals, ylim = c(0, 100),width = 1,space = c(1.1, rep(3, times = length(rev_vals) - 1)),col = 'magenta1',add = TRUE,yaxt = 'n')
    barplot(ctrl_vals, ylim = c(0, 100),width = 1,space = c(2.2, rep(3, times = length(ctrl_vals) - 1)),col = 'grey',add = TRUE,yaxt = 'n')
    axis(1, at = seq(0, (length(fwd_vals) - 1) * 4 + 1, by = 4)+1.5, labels = x_labels, cex.axis = 2)
    legend('topright',inset = c(-0.2, 0.1),legend = c(plot.nom,'Sim. Dist.'),col = c('blue', 'magenta1','grey'),fill = c('blue', 'magenta1','grey'),cex = 1.2,bty = "n")

    dev.off()

  }

  # time stamp for end of processing job
  show(paste0('Processing End:  ',Sys.time()))

}

}
