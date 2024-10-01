library(dplyr)

#Target.Data <- "/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/TARGETArrangednoxy.bed"
#target.file <- "/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/TARGETchrse.bed"
#Target.Data.df <- read.table(Target.Data, header = TRUE, sep = "\t", fill = TRUE)
#ref_qc <- GRanges(seqname = Target.Data.df$chr, IRanges(start = Target.Data.df$start, end = Target.Data.df$end))
#ref_qc <- Target.Data.df

ref_qc <- read.table('/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/NEWTAR.bed', 
                              header = TRUE, sep = "\t", fill = TRUE ) 
  
  
#..........................................

NEWAPP<- read.table('/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/Removalofzeros.bed', 
                    header = FALSE, sep = "\t", fill = TRUE )
#NORMALIZATION--Z-SCORE
reads<- NEWAPP$V5
expected_reads <- mean(reads)
z_normalized_ratio <- (reads - expected_reads)/sd(reads)
Z <- z_normalized_ratio

#write.csv(Z,paste0("/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/Viterbi-Code/my-vit-algorithm-using-for-thesis/Z-normalised1.bed"))

#aa <- median(Z)

# Define emission probabilities based on Z-normalized values
emission1 <- dnorm(Z,-0.05, 1)
emission1[emission1 == 0] <- min(emission1[emission1>0])
emission2 <- dnorm(Z, 0, 1)
emission2[emission2 == 0] <- min(emission2[emission2 > 0])
emission3 <- dnorm(Z, 0.05, 1)
emission3[emission3 == 0] <- min(emission3[emission3 > 0])

# Define transition probabilities

#My transition Probabilities I am using
#transition_probs <- matrix(c(0.9, 0.05, 0.05,
#                            0.05, 0.9, 0.05,
#                            0.05, 0.05, 0.9), nrow = 3, byrow = FALSE)
#.......................................................
transition_probs <- matrix(c(0.9, 0.05, 0.05,
                             0.05, 0.9, 0.05,
                             0.05, 0.05, 0.9), nrow = 3, byrow = FALSE)

                             

# Viterbi algorithm
n <- length(Z)
v <- matrix(NA, nrow = n, ncol = 3)
pointer <- matrix(NA, nrow = (n+1), ncol = 3)

# Initialization
v[1, 1] <- log(emission1[1]) + log(1 / 3)
v[1, 2] <- log(emission2[1]) + log(1 / 3)
v[1, 3] <- log(emission3[1]) + log(1 / 3)
#pointer[1, ] <- 0

# Recursion
for (i in 2:n) {
  for (j in 1:3) {
    max_prev <- max(v[i - 1, ] + log(transition_probs[, j]))
    v[i, j] <- log(get(paste0("emission", j))[i]) + max_prev
    pointer[i, j] <- which.max(v[i - 1, ] + log(transition_probs[, j]))
  }
}

# termination
max4=max(v[n,1],v[n,2],v[n,3])
if (v[n,1]==max4){
  pointer[n+1,1]=1
} else if (v[n,2]==max4) {
  pointer[n+1,1]=2
} else {
  pointer[n+1,1]=3
}

# traceback
traceback=pointer[n+1, 1]
result=rep(NA,n)
for (i in seq(n, 2, -1)){
  if (traceback==1){
    result[i-1]="deletion"
    traceback=pointer[i-1,1]
  } else if (traceback==2) {
    result[i-1]="neutral"
    traceback=pointer[i-1,2]
  } else {
    result[i-1]='duplication'
    traceback=pointer[i-1,3]
  }
}  

# Output result
cat("CNV calls:", result, "\n")

#Create Output-----------------------------------------------------------------
#finalcalls <- list()
for(chr in 1:24){
  delIndex = which(result=="deletion")
  dupIndex = which(result=="duplication")
    if( ( length( delIndex ) + length( dupIndex ) ) == 0 ) {
      finalcalls = NULL
    } else {
      stDel = delIndex[ !delIndex %in% ( delIndex + 1 ) ]
      edDel = delIndex[ !delIndex %in% ( delIndex - 1 ) ]
      stDup = dupIndex[ !dupIndex %in% ( dupIndex + 1 ) ]
      edDup = dupIndex[ !dupIndex %in% ( dupIndex - 1 ) ]
      st_exon = c( stDel, stDup )
      ed_exon = c( edDel, edDup )
      cnv_type = c( rep( "del", length( stDel ) ),
                rep( "dup", length( stDup ) ) )
      st_bp = ref_qc$start[ st_exon ]
      chr = ref_qc$chr[st_exon]
      ed_bp = ref_qc$end[ ed_exon ] # + ref_qc$width[ ed_exon ] - 1
      length_kb = ( ed_bp - st_bp) / 1000
       # if ( mode == "integer" ) copy_no = round( copy_no )
      finalcalls = data.frame( cnv_type, st_bp, ed_bp, chr,
                                  length_kb)
      colnames( finalcalls ) = c( 'cnv','st_bp','ed_bp', 'chr',
                                    'length_kb')
      }
}


write.table(finalcalls,paste0("/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/Viterbi-Code/my-vit-algorithm-using-for-thesis/resultsFIN.csv"))
