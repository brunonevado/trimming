countmissing<-function(ifas,iseq,istart,iend){
  x<-as.character(ifas[iseq,istart:iend])
  length(c(grep("-",x),grep("n",x)))
}

calcmeandist<-function(ifas,iseq,istart,iend,verbose=FALSE){
  if(verbose ==TRUE){
    message("calcmeandist,", iseq, ":",istart,"->",iend)
  }
  idiffs<-vector()
  for(i in 1:dim(ifas)[1]){
    if(i != iseq){
      idiffs<-c(idiffs,dist.dna(model="raw",ifas[c(iseq,i),istart:iend]))
    }
  }
  x<-mean(idiffs,na.rm=TRUE)
  if(is.na(x)){
    x<-0
  }
  x
}

trimfasta<-function(fas,thr=10,winlen=9,plots=FALSE,verbose=FALSE){

  # assumes input are CDS
  # will cut beginning and end of each sequence if too divergent from remaining sequences (> thr x average dist)
  # works in windows of length winlen, ignores missing data
  
  # 1. calculate mean distance of each sequence to all the remaining sequences
  diffs<-vector("numeric")
  meandiffs<-vector("numeric",dim(fas)[1])
  
  for(seq2check in 1:dim(fas)[1]){
    for(i in 1:dim(fas)[1]){
      if (i == seq2check){
       next
      }
    diffs<-c(diffs,dist.dna(model="raw",fas[c(seq2check,i),]))
    }
  meandiffs[seq2check]<-mean(diffs,na.rm=TRUE)
  diffs<-vector("numeric")
  }
  
  # this matrix will contain 1s for windows to remove
  distmat<-matrix(nrow=dim(fas)[1],ncol=ceiling(dim(fas)[2]/winlen))

  # 2. find front windows to remove
  for(seq2cut in 1:dim(fas)[1]){
    winn<-0
    for(wstart in seq(1,dim(fas)[2],winlen)){
      if(wstart + winlen > dim(fas)[2]){
        wend<-dim(fas)[2]
      }
      else{
        wend<-wstart + winlen -1
      }
      winn<-winn+1
      # skip if no data in window for seq of interest
      if(countmissing(fas,seq2cut,wstart,wend) == winlen){
        next
      }
      # if some data check divergence, if too high mask and repeat, otherwise finish
      if(verbose == TRUE) {
        message(seq2cut,":",wstart,"->",wend)
        }
      if(calcmeandist(fas,seq2cut,wstart,wend) > (thr*meandiffs[seq2cut])){
        if(verbose == TRUE) {
          message("pass")
          }
        
        if(plots==TRUE){
          par(mfrow=c(2,1))
          image(fas)
          abline(v=wstart,lwd=6,col="purple",lty=2)
          abline(v=wend,lwd=6,col="purple",lty=2)
          title(sub=paste("Masking sequence",rownames(fas)[seq2cut],"\nat positions:",wstart,"->",wend))
          image(fas[,wstart:wend])
        }
        distmat[seq2cut,winn]<-1
        next
      }
      else{
        break
      }
    }
  }

  # 3. find end windows to remove
  for(seq2cut in 1:dim(fas)[1]){
    winn<-dim(distmat)[2]+1
    for(wstart in rev(seq(1,dim(fas)[2],winlen))){
      if(wstart + winlen > dim(fas)[2]){
        wend<-dim(fas)[2]
      }
      else{
        wend<-wstart + winlen -1
      }
      winn<-winn-1
      # skip if no data in window for seq of interest
      if(countmissing(fas,seq2cut,wstart,wend) == winlen){
        next
      }
      # if some data check divergence, if too high mask and repeat, otherwise finish
      if(calcmeandist(fas,seq2cut,wstart,wend) > (thr*meandiffs[seq2cut])){
        if(plots==TRUE){
          par(mfrow=c(2,1))
          image(fas)
          abline(v=wstart,lwd=6,col="purple",lty=2)
          abline(v=wend,lwd=6,col="purple",lty=2)
          title(sub=paste("Masking sequence",rownames(fas)[seq2cut],"\nat positions:",wstart,"->",wend))
          image(fas[,wstart:wend])
        }
        distmat[seq2cut,winn]<-1
        next
      }
      else{
        break
      }
    }
  }
  
  # 4. remove windows
  for(seq2cut in 1:dim(fas)[1]){
    winn<-0
    for(wstart in seq(1,dim(fas)[2],winlen)){
      if(wstart + winlen > dim(fas)[2]){
        wend<-dim(fas)[2]
      }
      else{
        wend<-wstart + winlen -1
      }
      winn<-winn+1
      if(!is.na(distmat[seq2cut,winn])){
        fas<-as.character(fas)
        fas[seq2cut,wstart:wend]<-rep('n',wend-wstart+1)
        fas<-as.DNAbin(fas)
      }
    }
  }
  list(fas, length(which(distmat == 1)))
  
}

trimallfiles<-function(fastalist,thr=10,winlen=9,plots=FALSE,outsuffix){
  require(ape)
  fastafiles = readLines(fastalist)
  message(paste("<TrimAlign.R>: Read", length(fastafiles), "fasta files, threshold (x above average div) = ", thr,", window length =", winlen))
  nwindowsmasked<-0
  nfileschecked<-0
  for(i in 1:length(fastafiles)){
    if(!file.exists(fastafiles[i])){
      message("Cannot find input file ", fastafiles[i], ", skipping" )
      next
    }
    
    fas<-read.dna(fastafiles[i],format="f")
    x<-trimfasta(fas,thr=thr,winlen=winlen,plots=plots)
    nwindowsmasked<-nwindowsmasked+x[[2]]
    outfas<-paste(fastafiles[i],outsuffix,"fas",sep=".")
    write.dna(x[[1]],outfas,format="fasta",nbcol=-1,colsep="")
    nfileschecked<-nfileschecked+1
  }   
  message("Checked ", nfileschecked, " files, masked ", nwindowsmasked, " windows")
}


## script call

args <- commandArgs(TRUE)

if(length(args)!=4){
    stop( paste( "ERROR: wrong number of arguments supplied, expected 4 ( fastalist, threshold, winlen, outprefix ) but received ", length(args) ))
}

fastalist<-args[1]
thr<-as.numeric(args[2])
winlen<-as.numeric(args[3])
outsuffix<-args[4]

plots<-paste(sep="",outsuffix, ".trimmedEnds.pdf")
pdf(plots)
trimallfiles(fastalist=fastalist,thr=thr,outsuffix=outsuffix,winlen=winlen,plots=FALSE)
dev.off()


