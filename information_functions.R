
# return conditional entropy

entropy <- function(mat,norm=T){
  p.Sj=colSums(mat)/sum(colSums(mat))
  Hc = sum(sapply(1:ncol(mat), function(j){
    H.j = -sum(sapply(1:nrow(mat),function(i){
      p=mat[i,j]/sum(mat[,j]) #probability of species i, given signal j
      (p*log2(p)) #conditional entropy of species i, given signal j
    }),na.rm=TRUE) # sum conditional entropies of all species given signal j
    if (norm==T){
      H.j=H.j/log2(nrow(mat))}
    p.Sj[j]*H.j
  })# sum across all objects to get the conditional entropy Hc
  )
  return(Hc)
}


# return specific information

specific_information <- function(mat,norm=T){
  n=sum(colSums(mat))
  H <- sapply(1:nrow(mat), function(i){
    H.i = sum(sapply(1:ncol(mat),function(j){
      p.s <- sum(mat[,j])/n # prior probability of volatile Sj
      p.so <- mat[i,j]/sum(mat[i,]) # probability of signal Sj, given object o
      p.so*log2(p.so/p.s) # cf equation 3
    }),na.rm=T)
    if (norm==T){
      H.i <- H.i/log2(nrow(mat))}
  })
  return(H)
}
