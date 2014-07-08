####################################################
              ##REVERSE COMPLEMENT##
####################################################
myseq <- c("ATGCATTGGACGTTAG")
revComp <- function(x=myseq, mode="RC"){
  myvect <- NULL
  for(i in seq(along=(x[]))){
    if(mode=="C") {
      myvect[i] <- chartr("ATGC", "TACG", x[i])
    }
    if(mode=="RC" | mode=="R"){
      vecstring <- strsplit(x[i],"")[[1]]
      revv <- rev(vecstring)
      collapse <- paste(revv, collapse="")      
    if(mode=="RC")
      myvect <- chartr("ATGC", "TACG", collapse)
    }
  }
  return(myvect)
}
revComp(myseq)
####################################################
                ##FROM DATA SET##
####################################################

mydf <- DNAseq
myvec <- NULL
reverselist <- function(x, mode="RC"){
  for ( i in seq(along=mydf[,1])){
    myvec <- c(myvec,revComp(as.character(mydf[i,2])))   
  }
  return(myvec)
}

####################################################
              ##TRANSLATE TO PROTEIN##
####################################################
AAdf <- read.table(file="AA.txt", header=T, sep="\t")
AAv <- AAdf[,2]; names(AAv) <- AAdf[,1]
translate <- function(x=sequence, orf, list=AAv){
  sep <- function(x){
    seq <- gsub("(...)", "\\1_", x)
    seq <- unlist(strsplit(seq, "_"))
    seq <- seq[grep("^...$", seq)]
    return(seq)
  }
  if (orf==1){
    y <- sep(x)    
    trans=""
    for (i in 1:length(y)){
      seq_index <- y[i]
      trans <- paste(trans, as.character(AAv[seq_index]),sep="") 
      #print(trans)
    }
    #print(trans)
  }
  if (orf==2){
    y <- sep(x)
    x <- gsub("^.", "", x)
    trans=""
    for (i in 1:length(y)){
      seq_index <- y[i]
      trans <- paste(trans, as.character(AAv[seq_index]), sep="")
      }
    }
  if (orf==3){
    x <- gsub("^..", "", x)
    y <- sep(x)
    trans=""
    for (i in 1:length(y)){
      seq_index <- y[i]
      trans <- paste(trans, as.character(AAv[seq_index]), sep="")
    }
  }
return(trans)
}
sequence=c("TCACTACTTCGCTAA")
translate(sequence, 1, list=AAv)

