# ======================
# INTERNAL FUNCTIONS
# ======================

sumtrees<-function(intree,what="TL"){
  # given a tree, returns the total lengths or the tip or edges length (possible 
  # standardised by dividing by tree length)
  if(what == "TL"){
      # get total length
      return ( sum(intree$edge.length) )
  }
  else if(what == "in"){
      # get internal branch lengths
      return ( intree$edge.length[intree$edge[,2] > Ntip(intree)] )  
  }
  else if(what == "insd"){
      # get standardised internal branch lengths
      return ( intree$edge.length[intree$edge[,2] > Ntip(intree)] / sum(intree$edge.length) )  
  } 
  else if(what == "out"){
      # get tip lengths
      return ( intree$edge.length[intree$edge[,2] <= Ntip(intree)] )  
  }
  else if(what == "outsd"){
      # get standardised tip lengths
      return ( intree$edge.length[intree$edge[,2] <= Ntip(intree)] / sum(intree$edge.length) )  
  }   
  else{
      stop("ERROR: dont know what to get")
  }
}


sumtipsSD<-function(intree,ichars=20){
  # given a tree, returns the standardised (dividing by tree length) of each tip
  # assumes the first n chars in each sequence denote the individual (are common across trees)
  x<-list()
  for(i in 1:dim(intree$edge)[1]){
    if( intree$edge[i,2] <= Ntip(intree) ){
      tipname<-substr(intree$tip.label[ intree$edge[i,2]  ],1,ichars)
      tipsdlen<- intree$edge.length[i] / sum(intree$edge.length)
      x[[tipname]] <-tipsdlen
    }
  }
  # get standardised tip lengths
  return(x)
}

 
 
trimtips<-function(intree, thrlist,ichar=20,verbose=FALSE){
  # removes tips that are longer than standardised thr (i.e. length/TL > thr)
  # thr is a list, with a different thr for each species
  utr<-unroot(intree)
  tl<-sum(utr$edge.length)
  # scaled edge lengths
  edgelens<-utr$edge.length/tl
  # edges longer than tip threshold (only the right-hand side of edge)
  # has to be checked for each tip
  toremove<-vector("character")
  for(iedge in 1:dim(intree$edge)[1]){
    if( intree$edge[iedge,2] <= Ntip(intree)){
      tip_name<-intree$tip.label[intree$edge[iedge,2]]
      sp_name<-substr(tip_name,1,ichar)
      if( edgelens[iedge] > thrlist[[sp_name]] ){
        toremove<-c(toremove,tip_name)
      }
    }
  }
  

  # remove tips, but first checkthat at least 3 tips remain
  if( Ntip(utr) <= (length(toremove)+2) ){
    if (verbose == TRUE){
      message("Skipping tree: not enough tips after trimming tips")
      }
    return (NULL)
  }
  if(verbose==TRUE){
    message(paste("Removing" , length(rmtips), "tip(s)"))
  }
  utr2<-drop.tip(utr,toremove)

  return(utr2)
  }
  
  
    
hastoolongedge<-function(intree,thr,verbose=FALSE){
  # returns true if tree has any too long internal edge (length/TL > thr)
  if(verbose==TRUE){
    message(paste("Checking edges, tree:", write.tree(intree)))
  }
  utr<-unroot(intree)
  tl<-sum(utr$edge.length)
  
  # scaled edge lengths
  xedgelens<-utr$edge.length/tl
  # edges longer than tip threshold (only the right-hand side of edge)
  if(length(which(xedgelens> thr)) == 0){
    return (FALSE)
  }
  longedges<-utr$edge[which(xedgelens> thr),2]
  rmedges<-longedges [which(longedges > Ntip(utr))]
  if( length(rmedges> 0) ){
    return (TRUE)
  }
  else{
    return(FALSE)
  }
}

trimedges<-function(intree, thr,verbose=FALSE){
  # breaks tree at too long edges (length/TL > thr) and returns tree with more tips
  # works recursively until no too long edges remain
 if(verbose==TRUE){
   message(paste("trimming edges, tree:", write.tree(intree)))
   }

  utr<-unroot(intree)
  while(hastoolongedge(utr,thr) == TRUE){
    tl<-sum(utr$edge.length)
    # scaled edge lengths
    edgelens<-utr$edge.length/tl
    # edges longer than tip threshold (only the right-hand side of edge)
    longedges<-utr$edge[which(edgelens> thr),2]
    rmedges<-longedges [which(longedges > Ntip(utr))]
	str1<-extract.clade(utr,rmedges[1])
   	str2<-drop.tip(utr,str1$tip.label)
    # take tree with more tips (or random if same)
    if( Ntip(str1) > Ntip(str2) ){
      utr<-str1
    }
    else{
      utr<-str2
    }
    if(Ntip(utr) < 3){
      if(verbose){
      message("Skipping, not enough tips surviving")
      }
      return (NULL)
    }
    utr<-unroot(utr)
    if(verbose==TRUE){
        message(paste("Surviving subtree:",write.tree(utr)))
    }
  }
  return(utr)
  }

getTips<-function(intree){
  return(intree$tip.label)
}


# ======================
# MAIN
# ======================

trimtrees<-function(treelist,fastalist,thr=.99,explore=0,verbose=FALSE,outprefix="trimmed1", minSpecies=6, survival=FALSE){
  require(ape)
  treefiles = readLines(treelist)
  fastafiles = readLines(fastalist)
  if(length(fastafiles) != length(treefiles)){
    stop("Length of tree and fasta files list does not match")
  }
  n=length(fastafiles)
  trees<-vector("list", n)
  for(i in 1:n){
    trees[[i]]<-read.tree(treefiles[i])
  }
  class(trees)<-"multiPhylo"
  geneindx<-1:length(trees)
  
  # calculate thresholds
  utrs<-lapply(trees, unroot)
  tls<-lapply(utrs, sumtrees)
  inlens<-lapply(utrs, sumtrees,"in")
  inlenssd<-lapply(utrs, sumtrees,"insd")
  outlens<-lapply(utrs, sumtrees,"out")

  # getting SD tip length for each species
  x<-list()
  for(i in 1:length(utrs)){
    x<-c(x,sumtipsSD(utrs[[i]]))
    }
  # get species names (should be constant between trees)
  spnames<-names(table(as.factor(names(x))))
  sdthrsp<-list()
  for(i in 1:length(spnames)){
    sdthrsp[[spnames[i]]]<-quantile(unlist(x[grep(spnames[i],names(x))]),thr)
  }


  outlenssd<-lapply(utrs, sumtrees,"outsd")
  throut<-quantile(unlist(outlenssd),thr)
  thrin<-quantile(unlist(inlenssd),thr)
  # info
  message(paste("<TrimTrees.R>: Read", length(trees), "trees, threshold (internal edges) =", thrin))
  message(paste("Threshold (tips) = ", names(sdthrsp) ))
  message(paste("Threshold (tips) = ", sdthrsp ))
  
  # trim trees of too long tips
  #utrs_notips<-lapply(utrs, trimtips, throut)
  utrs_notips<-lapply(utrs, trimtips, sdthrsp)
  # summary and get rid of trees without enough branches
  exc1<-lapply(utrs_notips,is.null)
  Nexc_trees1<-length(which(exc1==TRUE))
  geneindx<-geneindx[!unlist(exc1)]
  utrs_pass1<-utrs_notips[!unlist(exc1)]
  total_tips<-sum(unlist(lapply(utrs,Ntip)))
  Nexc_tips1<-sum(unlist(lapply(utrs,Ntip)))-sum(unlist(lapply(utrs_pass1,Ntip)))
  message(paste("Excluded" , Nexc_trees1, "of",length(utrs)," trees due to trimming of tips"))
  message(paste("Excluded" , Nexc_tips1, "of",total_tips," tips due to trimming of tips"))

  # write excluded sequence names - to file, TODO
  alltips<-unlist(lapply(utrs,getTips))
  survivingtips1<-unlist(lapply(utrs_pass1,getTips))
  trimmedtips<-setdiff(alltips,survivingtips1)
  if(!is.null(outprefix)){
    out1<-paste(outprefix,"exctips.tips.txt",sep=".")
    write.table(trimmedtips,out1, quote=FALSE,row.names=FALSE, col.names=FALSE)
  }

  # break trees at too long internal edges
  utrs_inedges<-lapply(utrs_pass1, trimedges, thrin)

  # summary and get rid of trees without enough branches
  exc2<-lapply(utrs_inedges,is.null)
  Nexc_trees2<-length(which(exc2==TRUE))
  utrs_pass2<-utrs_inedges[!unlist(exc2)]
  geneindx<-geneindx[!unlist(exc2)]
  Nexc_tips2<- total_tips- Nexc_tips1-sum(unlist(lapply(utrs_pass2,Ntip)))
  message(paste("Excluded" , Nexc_trees2, "of",length(utrs_pass1)," trees due to trimming of internal edges"))
  message(paste("Excluded" , Nexc_tips2, "of",total_tips-Nexc_tips1," tips due to trimming of internal edges"))

  # write excluded sequence names to file
  alltips2<-unlist(lapply(utrs_pass1,getTips))
  survivingtips2<-unlist(lapply(utrs_pass2,getTips))
  trimmedtips2<-setdiff(survivingtips1,survivingtips2)

  if(!is.null(outprefix)){
    out2<-paste(outprefix,"exctips.edges.txt",sep=".")
    write.table(trimmedtips2,out2, quote=FALSE,row.names=FALSE, col.names=FALSE)
  }

  if( survival == TRUE){
    return (c(thr, length(utrs_pass2),length(survivingtips2),length(trees),total_tips))
  }
  
  # now read and trim alignments. show some examples if explore
  toofewsp<-0
  if( explore > 0 ){
    # this ensures we look only at genes that survived trimming 
    check<-sample(geneindx,explore)
    
    for(i in check){
      #message(i)
      original_tr<-trees[[i]]
      original_fas<-read.dna(fastafiles[i],format="f")
      surviving_tr<- utrs_pass2[[which(geneindx==i)]]
      surviving_fas<-original_fas[surviving_tr$tip.label,]
      
      par(mfrow=c(2,2))
      image(original_fas)
      plot(original_tr)
      add.scale.bar()
      image(surviving_fas)
      plot(surviving_tr,main=paste("N tips removed:", Ntip(original_tr) - Ntip(surviving_tr)))
      add.scale.bar()
     # readline(prompt="Press [enter] to continue")
    }
  }
  else{
      for(i in geneindx){
        surviving_tr<- utrs_pass2[[which(geneindx==i)]]
        if(Ntip(surviving_tr) < minSpecies){
          toofewsp<-toofewsp+1
          next
        }
        original_tr<-trees[[i]]
        original_fas<-read.dna(fastafiles[i],format="f")
        surviving_fas<-original_fas[surviving_tr$tip.label,]
        if (length( surviving_tr$tip.label %in% original_tr$tip.label ) != Ntip(surviving_tr)){
          stop("oops - mismatch between input and surviving trees' tip names!")
        }
        
        outtree<-paste(treefiles[i],outprefix,"tre",sep=".")
        outfas<-paste(fastafiles[i],outprefix,"fas",sep=".")
        write.tree(surviving_tr,outtree)
        write.dna(surviving_fas,outfas,format="fasta",nbcol=-1,colsep="")
        
      }
  message(paste("Number of trees excluded due to few sequences:", toofewsp))
  }
 
 
}



## script call

args <- commandArgs(TRUE)

if(length(args)!=6){
    stop( paste( "ERROR: wrong number of arguments supplied, expected 6 ( treelist, fastalist, threshold, outprefix, minSpecies, exploreTreeN ) but received ", length(args) ))
}

treelist<-args[1]
fastalist<-args[2]
thr<-as.numeric(args[3])
thr<-thr/100
outprefix<-args[4]
minSpecies<-as.numeric(args[5])
explore<-as.numeric(args[6])

if(explore > 0){
  plots<-paste(sep="",outprefix, ".trimmedTrees.pdf")
  pdf(plots)
}
trimtrees(treelist,fastalist,thr=thr,outprefix=outprefix, minSpecies=minSpecies, explore)
if(explore > 0){
  dev.off()
}
