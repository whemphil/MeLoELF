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
                    sam.file=list.files(path = mdir,pattern = '*.sam'),
                    crunch.too=T,
                    process=T,
                    pre.ligated=F,
                    methyl.type='B',
                    thresh.meth='BM',
                    BM.strand.data='r',
                    target_fwd_auc=0.9,
                    read.length=c(10,3000),
                    completeness=0.9,
                    matching=0.9,
                    align.file='align.RData',
                    processed.file='processed.RData',
                    seq.file='sequences.txt',
                    melo.file='methlocations.txt',
                    meca.file='methcalls.txt',
                    plot_title=''){


#######################################################################
### Script (don't edit past here)
#######################################################################

setwd(mdir)

#######################
## Load necessary non-base packages
#######################

library(stringr)
library(foreach)
library(doParallel)

#######################
## Create custom functions for later analysis
#######################

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
    if(sum(s1>1)==0){
      return(NULL)
    }
    s2=range(which(s1>1))
    results=rep(F,times=length(data))
    results[s2[1]:s2[2]]=T
    return(results)
  }

  test.0=str_extract_all(read,boundary("character"))[[1]] # takes the polymer sequence and converts it from a single string into a vector of 1 base per value
  lengths.of.reads=length(test.0)
  if((length(Cm)+length(Chm))!=length(C.key) | length(test.0)<min(read.length) | length(test.0)>max(read.length)){
    DATA='blank' # bypasses polymers outside desired length range
  } else {
    test.Chm=rep(NA,times=length(test.0)) # creates temporary  vector to store hydroxy methyl scores for polymer, with NA as default value
    test.Cm=test.Chm
    test.Chm[which(test.0=='C')] = 0 # sets the default methylation score for all Cs in the polymer to zero
    test.Cm[which(test.0=='C')] = 0
    try(test.Chm[which(test.0=='C')[Chm]]<-C.key[1:(0.5*length(C.key))]) # writes over methylation scores for annotated Cs in the polymer
    try(test.Cm[which(test.0=='C')[Cm]]<-C.key[(0.5*length(C.key)+1):length(C.key)])
    test.Chm=c(rep(NA,times=length(FWD)),test.Chm,rep(NA,times=length(FWD))) # adds empty flanks to the edges of methyl score vectors, in preparation for sliding alignment
    test.Cm=c(rep(NA,times=length(FWD)),test.Cm,rep(NA,times=length(FWD)))
    test.00=c(rep("x",times=length(FWD)),test.0,rep("x",times=length(FWD))) # adds empty flanks to edges of sequence vector, in preparation for sliding alignment
    refs=rep(0,times=length(test.00))
    test.1=test.00
    testy.Chm=test.Chm
    testy.Cm=test.Cm

    score.FWD=rep(0,times=length(test.1)+1-length(FWD)) # creates temporary vectors to hold the alignment scores for the FWD and REV sequences to the reference polymer
    score.REV=rep(0,times=length(test.1)+1-length(REV))

    NN=round(length(test.0)/length(FWD)+0.4)+2 # sets the maximum number of fragments allowed to be mapped to the polymer
    FWD.mat=matrix(FWD,nrow=NN,ncol = length(FWD),byrow = T)
    REV.mat=matrix(REV,nrow=NN,ncol = length(REV),byrow = T)
    FWD.align=matrix('x',nrow=NN,ncol = length(FWD)); colnames(FWD.align)<-paste0(FWD,'.',1:length(FWD))
    REV.align=matrix('x',nrow=NN,ncol = length(REV)); colnames(REV.align)<-paste0(REV,'.',1:length(REV))
    FWD.Chm=matrix(NA,nrow=NN,ncol = length(FWD)); colnames(FWD.Chm)<-paste0(FWD,'.',1:length(FWD))
    REV.Chm=matrix(NA,nrow=NN,ncol = length(REV)); colnames(REV.Chm)<-paste0(REV,'.',1:length(REV))
    FWD.Cm=matrix(NA,nrow=NN,ncol = length(FWD)); colnames(FWD.Cm)<-paste0(FWD,'.',1:length(FWD))
    REV.Cm=matrix(NA,nrow=NN,ncol = length(REV)); colnames(REV.Cm)<-paste0(REV,'.',1:length(REV))
    COUNTER=0

    FWD.refs=matrix(0,nrow=NN,ncol = length(FWD))
    REV.refs=matrix(0,nrow=NN,ncol = length(REV))

    # loop for sliding reference sequence along read to calculate positional alignment scores
    for (i in 1:length(score.FWD)){

      score.FWD[i]=mean(FWD==test.1[i:(i+length(FWD)-1)])-mean(FWD!=test.1[i:(i+length(FWD)-1)] & "x"!=test.1[i:(i+length(FWD)-1)])*0.1 # quantifies the alignment quality of FWD reference fragment at various positions along the polymer
      score.REV[i]=mean(REV==test.1[i:(i+length(REV)-1)])-mean(REV!=test.1[i:(i+length(REV)-1)] & "x"!=test.1[i:(i+length(REV)-1)])*0.1
      # The following sections prioritize the perfect reference sequence matches in the polymer, recording them as the first mapped fragments
      if (score.FWD[i]==1){
        COUNTER=COUNTER+1
        FWD.align[COUNTER,]=test.1[i:(i+length(FWD)-1)]
        FWD.Chm[COUNTER,]=test.Chm[i:(i+length(FWD)-1)]
        FWD.Cm[COUNTER,]=test.Cm[i:(i+length(FWD)-1)]
        FWD.refs[COUNTER,]=i:(i+length(FWD)-1)
        test.1[i:(i+length(FWD)-1)]=rep("x",times=length(FWD))
        test.Chm[i:(i+length(FWD)-1)]=rep(NA,times=length(FWD))
        test.Cm[i:(i+length(FWD)-1)]=rep(NA,times=length(FWD))
      }
      if (score.REV[i]==1){
        COUNTER=COUNTER+1
        REV.align[COUNTER,]=test.1[i:(i+length(REV)-1)]
        REV.Chm[COUNTER,]=test.Chm[i:(i+length(REV)-1)]
        REV.Cm[COUNTER,]=test.Cm[i:(i+length(REV)-1)]
        REV.refs[COUNTER,]=i:(i+length(REV)-1)
        test.1[i:(i+length(REV)-1)]=rep("x",times=length(REV))
        test.Chm[i:(i+length(REV)-1)]=rep(NA,times=length(REV))
        test.Cm[i:(i+length(REV)-1)]=rep(NA,times=length(REV))
      }

    }

    # loop for performing additional reference alignments with scoring, until sufficient read coverage is achieved
    for (j in (COUNTER+1):NN){

      # terminates further alignment attempts when less than 10 bases in the polymer remain unmapped
      if(sum(test.1!="x")<10){
        break
      }

      score.FWD=rep(0,times=length(test.1)+1-length(FWD))
      score.REV=rep(0,times=length(test.1)+1-length(FWD))

      # scoring loop
      for (i in 1:length(score.FWD)){
        score.FWD[i]=mean(FWD==test.1[i:(i+length(FWD)-1)])-mean(FWD!=test.1[i:(i+length(FWD)-1)] & "x"!=test.1[i:(i+length(FWD)-1)])*0.1
        score.REV[i]=mean(REV==test.1[i:(i+length(REV)-1)])-mean(REV!=test.1[i:(i+length(REV)-1)] & "x"!=test.1[i:(i+length(REV)-1)])*0.1
      }

      # Looks for the best-matching fragment alignments to the polymer, records that section as mapped, then repeats until the polymer is sufficiently mapped
      FWD.max=which.max(score.FWD)
      REV.max=which.max(score.REV)
      if(max(score.FWD)>=max(score.REV)){
        if(is.null(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))==T){
          break
        }
        FWD.align[j,which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))]=test.1[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]
        FWD.Chm[j,which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))]=test.Chm[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]
        FWD.Cm[j,which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))]=test.Cm[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]
        FWD.refs[j,which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))]=(FWD.max:(FWD.max+length(FWD)-1))[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))]
        test.Chm[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]=NA
        test.Cm[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]=NA
        test.1[which(bounds(test.1[FWD.max:(FWD.max+length(FWD)-1)]==FWD))+FWD.max-1]="x"
      }
      if(max(score.FWD)<max(score.REV)){
        if(is.null(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))==T){
          break
        }
        REV.align[j,which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))]=test.1[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]
        REV.Chm[j,which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))]=test.Chm[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]
        REV.Cm[j,which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))]=test.Cm[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]
        REV.refs[j,which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))]=(REV.max:(REV.max+length(REV)-1))[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))]
        test.Chm[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]=NA
        test.Cm[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]=NA
        test.1[which(bounds(test.1[REV.max:(REV.max+length(REV)-1)]==REV))+REV.max-1]="x"
      }

    }

    # for each reference sequence mapped to the polymer, get indices of the fragment within the polymer
    REV.refs[REV.align!=REV.mat & REV.align!="N"]=-1*REV.refs[REV.align!=REV.mat & REV.align!="N"]
    REV.refs[REV.align=="x"]=0
    FWD.refs[FWD.align!=FWD.mat & FWD.align!="N"]=-1*FWD.refs[FWD.align!=FWD.mat & FWD.align!="N"]
    FWD.refs[FWD.align=="x"]=0

    # calculate quality statistics for fragment alignments to read
    FEF.fwd=1-rowSums(FWD.align=='x')/length(FWD)
    FEF2.fwd=rowSums(FWD.align==FWD.mat)/rowSums(FWD.align!='x')
    FEF3.fwd=rowSums(FWD.align=='x')/length(FWD)
    FEF.rev=1-rowSums(REV.align=='x')/length(REV)
    FEF2.rev=rowSums(REV.align==REV.mat)/rowSums(REV.align!='x')
    FEF3.rev=rowSums(REV.align=='x')/length(REV)
    comp=rep(NA,times=NN)
    comp[FEF3.fwd!=1]=FEF.fwd[FEF3.fwd!=1]
    comp[FEF3.rev!=1]=FEF.rev[FEF3.rev!=1]
    match=rep(NA,times=NN)
    match[FEF3.fwd!=1]=FEF2.fwd[FEF3.fwd!=1]
    match[FEF3.rev!=1]=FEF2.rev[FEF3.rev!=1]
    id=rep(NA,times=NN)
    id[FEF3.fwd!=1]='FWD'
    id[FEF3.rev!=1]='REV'
    quality=cbind(comp,match,id)

    # save relevant data from read to subset of stable variable
    DATA=list('FWD.align'=FWD.align,'REV.align'=REV.align,'FWD.Chm'=FWD.Chm,'FWD.Cm'=FWD.Cm,'REV.Chm'=REV.Chm,'REV.Cm'=REV.Cm,'Q'=quality,'FWD.I'=abs(FWD.refs),'REV.I'=abs(REV.refs))

  }

  return(DATA)

}

# Mixed beta thresholding
BM.thresh <- function(data.actual.fwd,data.actual.rev,met='RSS',p=0.95,set=BM.strand.data){
  if(set=='b'){
    data=as.numeric(na.omit(c(data.actual.fwd[data.actual.fwd>0 & data.actual.fwd<1],data.actual.rev[data.actual.rev>0 & data.actual.rev<1])))
  }
  if(set=='f'){
    data=as.numeric(na.omit(c(data.actual.fwd[data.actual.fwd>0 & data.actual.fwd<1])))
  }
  if(set=='r'){
    data=as.numeric(na.omit(c(data.actual.rev[data.actual.rev>0 & data.actual.rev<1])))
  }
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
    cdf=data.frame(x=fit.dens$x,y=mean(diff(fit.dens$x))*cumsum(fit.dens$y))
    sep.point=fit.dens$x[which.min(fit.dens$y)]
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
  try(fit.betasSIN <- optim(par = init.par.est.single(data,fit.dens),fn = reg.beta.single,method = "L-BFGS-B"))
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
  dBIC=log(length(data))*length(fit.betas$par)+2*reg.beta(fit.betas$par,meth = 'NLL')
  sBIC=log(length(data))*length(fit.betasSIN$par)+2*reg.beta.single(par = fit.betasSIN$par,meth = 'NLL')
  dBeta.share=dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2])-dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]);dBeta.share[dBeta.share<0]=0;dBeta.share=(1-mean(diff(fit.dens$x),na.rm = T)*sum(dBeta.share))
  if(dBIC<sBIC & dBeta.share<0.5 & abs(0.5-fit.betas$par[5])<0.4){
    beta.means=as.numeric((fit.betas$par[c(1,3)])/(fit.betas$par[c(1,3)]+fit.betas$par[c(2,4)]))
    if(which.max(beta.means)==1){
      rel.lik=((fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]))/((1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4])))
    }
    if(which.max(beta.means)==2){
      rel.lik=1/((fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]))/((1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4])))
    }
    thresh=qbeta(p,fit.betas$par[2*(which.min(beta.means)-1)+1],fit.betas$par[2*(which.min(beta.means)-1)+2])
    thresh2=fit.dens$x[min(which(rel.lik>1 & fit.dens$x>min(beta.means)))]
    plot(fit.dens,ylim=c(0,max(fit.dens$y)),main='Beta Unmixing Threshold',xlab='Methyl Score');lines(fit.dens$x,fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2])+(1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]),col='red');abline(v=thresh,col='green',lwd=2,lty='dashed');abline(v=thresh2,col='blue',lwd=2,lty='dashed');legend('topright',legend = c('Data','Model','Threshold','Threshold 2'),col=c('black','red','green','blue'),fill=c('black','red','green','blue'))
    lines(fit.dens$x,fit.betas$par[5]*dbeta(fit.dens$x,shape1 = fit.betas$par[1],shape2 = fit.betas$par[2]),lty='dotted',col='purple',lwd=2);lines(fit.dens$x,(1-fit.betas$par[5])*dbeta(fit.dens$x,shape1 = fit.betas$par[3],shape2 = fit.betas$par[4]),lty='dotted',col='purple',lwd=2)
    #thresh=thresh2
  }
  if(sBIC<dBIC | dBeta.share>=0.5 | abs(0.5-fit.betas$par[5])>0.4){
    plot(fit.dens,ylim=c(0,max(fit.dens$y)),main='Beta Unmixing Threshold',xlab='Methyl Score');lines(fit.dens$x,dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]),col='red')
    show('WARNING! -- Methyl score values appear to primarily belong to a single distribution, and its corresponding methylation state is needed for thresholding -- attempting auto-assignment...')
    if(fit.dens$x[which.max(fit.dens$y)]>0.5){
      usr.input='p'
      show('...distribution inferred to correspond to methylated CpGs.')
    }else{
      if(fit.dens$x[which.max(fit.dens$y)]<0.15){
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
    }
    if(usr.input=='p'){
      neg.dens=fit.dens$y-dbeta(fit.dens$x,shape1 = fit.betasSIN$par[1],shape2 = fit.betasSIN$par[2]);neg.dens[neg.dens<0 | fit.dens$x>fit.dens$x[which.max(fit.dens$y)]]=0
      thresh=fit.dens$x[min(which((cumsum(neg.dens)/sum(neg.dens))>p))]
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
  bin_mat <- ifelse(mat_num < threshold, 0, 1)
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

# Save relevant parameters
PAR.list=ls()

#######################
## Read Parsing and Fragment Mapping
#######################


if(crunch.too){

  # warning for Windows PC users
  if(Sys.info()['sysname']=='Windows'){
    show('WARNING! -- The MeLoELF package was designed with Unix systems in mind, and may not work on Windows PCs.')
    show('The sam -> txt file processing invokes bash/awk functionality, so at minimum, the sequence, methlocations, and methcalls txt files must be supplied manually to the mdir directory.')
  }

  # time stamp for beginning of alignment job
  show(paste0('Align Start:  ',Sys.time()))

  # process relevant sam file information into individual txt files using bash/awk
  try(system(paste0("awk '{print $10}' ",getwd(),"/",sam.file," > ",getwd(),"/",seq.file)))
  try(system(paste0("awk '{print $28}' ",getwd(),"/",sam.file," > ",getwd(),"/",melo.file)))
  try(system(paste0("awk '{print $29}' ",getwd(),"/",sam.file," > ",getwd(),"/",meca.file)))
  #
  raw=read.csv(file = seq.file,header = F) # load sequences from pre-processed file
  raw.2=as.matrix(read.csv(melo.file,header = F,sep = ";")[,]) # load CpG indices from pre-processed file
  raw.3=read.csv(meca.file,header = F,sep = ";") # load methyl and hydroxy-methyl scores from pre-processed file

  # convert input reference sequences into vector of bases
  FWD=str_extract_all(parent,boundary("character"))[[1]]
  REV=str_extract_all(target,boundary("character"))[[1]]

  # load CpG indices, hydroxy-methyl scores, and methyl scores into alignable arrays
  Chm=str_split(raw.2[,1],',')
  Cm=str_split(raw.2[,2],',')
  C.key=str_split(raw.3[,1],',')
  for (i in 1:length(C.key)){
    Chm[[i]]=as.numeric(Chm[[i]][-1])
    Cm[[i]]=as.numeric(Cm[[i]][-1])
    C.key[[i]]=round(as.numeric(C.key[[i]][-1])/256,2)
    Chm[[i]]=cumsum(Chm[[i]]+1)
    Cm[[i]]=cumsum(Cm[[i]]+1)
  }

  # create empty variables for data storage and indexing
  READS=0
  lengths.of.reads=rep(NA,times=nrow(raw))

  # loop to perform parallelized alignments with scoring on a per-read basis
  registerDoParallel(detectCores())
  DATA <- foreach(seq = 1:nrow(raw)) %dopar% {

    map.fragments(read=raw[seq,1],Cm[[seq]],Chm[[seq]],C.key[[seq]],read.length=read.length,FWD=FWD,REV=REV)

  }
  # get lengths of reads
  lengths.of.reads <- foreach(seq = 1:nrow(raw),.combine = c) %dopar% {

    length(str_extract_all(raw[seq,1],boundary("character"))[[1]])

  }
  # get fragment numbers
  READS <- foreach(seq = 1:nrow(raw),.combine = c) %dopar% {

    round(length(str_extract_all(raw[seq,1],boundary("character"))[[1]])/length(FWD)+0.4)+2

  }

  # add final useful data to end of stable variable
  DATA[['N']]=sum(as.numeric(READS))
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

  # generate empty matrices to consolidate data
  FWD.Chm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.Chm)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.Chm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.Chm)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  FWD.Cm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.Cm)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.Cm=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.Cm)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  Q.reads=matrix(NA,nrow=DATA[['N']],ncol = 4); colnames(Q.reads)<-c('comp','match','id','polyN')
  FWD.index=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['FWD']])); colnames(FWD.index)<-paste0(DATA[['FWD']],'.',1:length(DATA[['FWD']]))
  REV.index=matrix(NA,nrow=DATA[['N']],ncol = length(DATA[['REV']])); colnames(REV.index)<-paste0(DATA[['REV']],'.',1:length(DATA[['REV']]))
  remapped.reads=as.data.frame(matrix(NA,nrow=length(DATA[['RSs']]),ncol=4)); colnames(remapped.reads) <- c('Length','Strand','Mcomp','Mmatch')
  COUNTER=1

  # loop for pulling alignment data out of large DATA containers and consolidating it
  for (i in 1:(length(DATA)-5)){

    if(DATA[i]=='blank'){
      next
    }

    FWD.Chm[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=DATA[[i]][['FWD.Chm']]
    REV.Chm[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=DATA[[i]][['REV.Chm']]
    FWD.Cm[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=DATA[[i]][['FWD.Cm']]
    REV.Cm[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=DATA[[i]][['REV.Cm']]
    Q.reads[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),1:3]=DATA[[i]][['Q']]
    Q.reads[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),4]=i
    if(pre.ligated){
      FWD.index[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=abs(DATA[[i]][['FWD.I']])-length(DATA[['FWD']])
      REV.index[COUNTER:(COUNTER+nrow(DATA[[i]][['Q']])-1),]=abs(DATA[[i]][['REV.I']])-length(DATA[['REV']])
      remapped.reads$Length[i]=DATA[['RLs']][i]
      remapped.reads$Mcomp[i]=mean(as.numeric(DATA[[i]]$Q[,1]),na.rm = T)
      remapped.reads$Mmatch[i]=mean(as.numeric(DATA[[i]]$Q[,2]),na.rm = T)
      if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
        remapped.reads$Strand[i]='FWD'
      }
      if(sum(DATA[[i]]$Q[,3]=='REV',na.rm = T)==sum(!is.na(DATA[[i]]$Q[,3]))){
        remapped.reads$Strand[i]='REV'
      }
      if(sum(DATA[[i]]$Q[,3]=='FWD',na.rm = T) < sum(!is.na(DATA[[i]]$Q[,3])) & sum(DATA[[i]]$Q[,3]=='REV',na.rm = T) < sum(!is.na(DATA[[i]]$Q[,3]))){
        remapped.reads$Strand[i]=sum(as.numeric(DATA[[i]]$Q[which(DATA[[i]]$Q[,3]=='REV'),1]),na.rm = T)/sum(as.numeric(DATA[[i]]$Q[,1]),na.rm = T)
      }
    }
    COUNTER=COUNTER+nrow(DATA[[i]]$Q)

  }
  if(pre.ligated){
    FWD.index[FWD.index<0]=NA
    REV.index[REV.index<0]=NA
  }

  if(pre.ligated){

    # compile data for REV polymers with sufficient mapping quality
    polymer.ids=which(remapped.reads$Strand=='REV' & remapped.reads$Mcomp>=completeness & remapped.reads$Mmatch>=matching)
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
    polymer.ids.fwd=which(remapped.reads$Strand=='FWD' & remapped.reads$Mcomp>=completeness & remapped.reads$Mmatch>=matching)
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

  #clean up data sets
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

    fwd.binary <- binarize_matrix(polymer.actual.fwd, thresh)
    rev.binary <- binarize_matrix(polymer.actual.rev, thresh)
  }
  if(!pre.ligated){
    # Get threshold values
    if(thresh.meth=='AUC'){
      thresh <- find_thresh_for_auc(data.actual.fwd, data.actual.rev)
    }
    if(thresh.meth=='BM'){
      thresh=BM.thresh(data.actual.fwd,data.actual.rev)
    }

    fwd.binary <- binarize_matrix(data.actual.fwd, thresh)
    rev.binary <- binarize_matrix(data.actual.rev[,ncol(data.actual.rev):1], thresh)
  }
  ctrl.binary=rev.binary; ctrl.binary[!is.na(ctrl.binary)]=0
  ctrl.binary[sample(which(!is.na(rev.binary)),sum(rev.binary,na.rm = T),prob = rep(colSums(rev.binary,na.rm = T)/colSums(!is.na(rev.binary)),times=colSums(!is.na(rev.binary))))]=1

  fwd.binary.surv <- get_survival_data(fwd.binary)
  rev.binary.surv <- get_survival_data(rev.binary)
  ctrl.binary.surv <- get_survival_data(ctrl.binary)

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

  save(matrices,thresh, file="methylation_matrices.RDS")
  write.csv(summary, "survival_summary.csv")

  # quality control analyses
  if(!pre.ligated){
    read.filt=nrow(data.actual.rev)/nrow(REV.Cm[rowSums(!is.na(REV.Cm))>0,])
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
        f53m.rev[i]=min(which(rev.binary[i,]==1))
        f35m.rev[i]=max(which(rev.binary[i,]==1))
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
    png('QCgraphs.png', height = 2650, width = 2800, res=300)
    par(mfrow=c(2,1),mar=c(3,6,3,1))
    #
    FWD.Cm.pdfs=list(NULL)
    REV.Cm.pdfs=list(NULL)
    for (i in 1:length(FWD.sites)){
      FWD.Cm.pdfs[[i]]=density(data.actual.fwd,na.rm = T,from = 0,to = 1)
      REV.Cm.pdfs[[i]]=density(data.actual.rev,na.rm = T,from = 0,to = 1)
    }
    plot(NULL,NULL,ylim=c(0,1),xlim=c(0,length(FWD.sites)),main = 'Methyl Score Distributions',cex.axis = 1.6,ylab = 'Methyl Score',cex.lab=2,cex.main=2,xaxt='n',xlab='')
    for(i in 1:length(FWD.sites)){
      lines(i-0.725-FWD.Cm.pdfs[[i]][['y']]/max(FWD.Cm.pdfs[[i]][['y']])*0.2,FWD.Cm.pdfs[[i]][['x']],col='blue')
      lines(i-0.725+FWD.Cm.pdfs[[i]][['y']]/max(FWD.Cm.pdfs[[i]][['y']])*0.2,FWD.Cm.pdfs[[i]][['x']],col='blue')
      lines(i-0.275-REV.Cm.pdfs[[i]][['y']]/max(REV.Cm.pdfs[[i]][['y']])*0.2,REV.Cm.pdfs[[i]][['x']],col='red')
      lines(i-0.275+REV.Cm.pdfs[[i]][['y']]/max(REV.Cm.pdfs[[i]][['y']])*0.2,REV.Cm.pdfs[[i]][['x']],col='red')
    }
    points(1:length(FWD.sites)-0.725,apply(data.actual.fwd,2,median,na.rm=T),col='purple',pch='-',cex=5)
    points(1:length(REV.sites)-0.275,apply(data.actual.rev,2,median,na.rm=T),col='purple',pch='-',cex=5)
    axis(side = 1,at = c(1:length(FWD.sites))-0.5,labels = paste0(DATA[['FWD']][FWD.sites],'p',DATA[['FWD']][FWD.sites+1],'-',FWD.sites+0.5),cex.axis=1.6)
    abline(h=thresh,lty='solid',col='grey',lwd=2)
    legend('topright',legend = c('FWD','REV'),col = c('blue','red'),fill = c('blue','red'),cex=1.2)
    #
    plot(NULL,NULL,xlim=c(0,5*length(FWD.sites)+2),ylim=c(0,1),xaxt='n',ylab='Fraction of Product',main='Methyl Location Distributions',cex.main=2,cex.lab=1.5,cex.axis=1.5,xlab='')
    axis(side = 1,at = c(0:6)*5+2,labels = paste0(DATA[['FWD']][FWD.sites],'p',DATA[['FWD']][FWD.sites+1],'-',FWD.sites+0.5),cex.axis=1.1)
    abline(h=fMs.rev,col='red',lty='dashed',lwd=2)
    abline(h=fSs.rev,col='purple',lty='dashed',lwd=2)
    abline(h=read.filt,col='orange',lty='dotted',lwd=2)
    points(c(0:6)*5+0.5,fCpGs.fwd,type='h',pch=22,lwd=15,col=1)
    points(c(0:6)*5+1.5,rev(fCpGs.rev),type='h',pch=22,lwd=15,col=2)
    points(c(0:6)*5+2.5,rev(f53m.rev),type='h',pch=22,lwd=15,col=3)
    points(c(0:6)*5+3.5,rev(f35m.rev),type='h',pch=22,lwd=15,col=4)
    legend('topright',legend=c('SYNTH Methyls','CAT Methyls',"5'->3' Start","3'->5' Start"),col=1:4,fill=1:4,cex=0.7)
    text(x=5*length(FWD.sites)-1,y=c(0.72,0.62,0.52),pos = 4,col = c('red','purple','orange'),cex = 1.0,labels = c('% CpG','% Sub.','% Frag.'))
    dev.off()
  }

  wd <- getwd()

  # plotting survival analysis
  png(paste0(wd, "/methylation_survival.png"), height = 1320, width = 2800, res=300)
  par(mfrow=c(1,1),mar=c(5,6,3,8))
  plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,100),
       xlab="Probability of successive methylation", ylab="Percentage of reads",
       main=paste(plot_title, " - Methylation survival"), cex.main = 2, cex.axis=2, cex.lab = 2)
  #main=expression(paste(Delta, "351 - Methylation survival")), cex.axis=1.2, cex.lab = 1.3)
  lines(fwd.surv.plot$x, fwd.surv.plot$y, type="s", col="blue", lwd=2)
  lines(rev.surv.plot$x, rev.surv.plot$y, type="s", col="magenta1", lwd=2)
  lines(ctrl.surv.plot$x, ctrl.surv.plot$y, type="s", col="grey", lwd=2)
  legend('topright', legend = c(paste0('SYNTH (AUC=',round(fwd.auc, digits=2),')'), paste0('CAT (AUC=',round(rev.auc, digits=2),')'), paste0('Sim. Dist. (AUC=',round(ctrl.auc, digits=2),')')),col = c('blue','magenta1','grey'),fill = c('blue','magenta1','grey'),cex=1.2, bty="n")
  #grid()
  #dev.copy2pdf(file="methylation_survival.pdf", height = 5, width = 7)
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

    png(paste0(wd, "/methyl_distribution.png"), height = 1320, width = 2800, res=300)
    par(mfrow=c(1,1), mar=c(5,6,3,8), xpd=TRUE)


    suppressWarnings(
      barplot(fwd_vals,ylim = c(0, 100),width = 1,space = c(0, rep(3, times = length(fwd_vals) - 1)),col = 'blue',cex.main = 2,cex.axis = 2,ylab = 'Percent of reads',xlab = 'Number of 5mCs',cex.lab = 2,main = paste(plot_title, " - Methylation distribution")))
    barplot(rev_vals, ylim = c(0, 100),width = 1,space = c(1.1, rep(3, times = length(rev_vals) - 1)),col = 'magenta1',add = TRUE,yaxt = 'n')
    barplot(ctrl_vals, ylim = c(0, 100),width = 1,space = c(2.2, rep(3, times = length(ctrl_vals) - 1)),col = 'grey',add = TRUE,yaxt = 'n')
    axis(1, at = seq(0, (length(fwd_vals) - 1) * 4 + 1, by = 4)+1.5, labels = x_labels, cex.axis = 2)
    legend('topright',inset = c(-0.2, 0.1),legend = c('SYNTH', 'CAT','Sim. Dist.'),col = c('blue', 'magenta1','grey'),fill = c('blue', 'magenta1','grey'),cex = 1.2,bty = "n")

    dev.off()

  }

  # time stamp for end of processing job
  show(paste0('Processing End:  ',Sys.time()))

}

}
