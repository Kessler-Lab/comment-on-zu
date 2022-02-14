### entropy functions from Zu et al 2020:

## Calculate entropy: 
# P_vector is a probability vector 
entropy <- function(P_vector){
  entropy <- 0
  for (i in 1:length(P_vector)) {
    if(P_vector[i] != 0) {  #zero check in class
      entropy0 <- -(P_vector[i]*log(P_vector[i], base = length(P_vector))) # for binary data set 
    }else{
      entropy0 <- 0
    }
    entropy <- entropy + entropy0
  }
  return(entropy)  #dont forget return!
}

# for average individual VOC when feeding a matrix
H_A_VOC <- function(x1) {
  x <- x1[, colSums(x1)!=0] # remove columns that only contains 0
  NI <- nrow(x) # number of species
  NJ <- ncol(x) # number of VOCs
  if (is.null(NJ)) {
    Hn_V <- 1
    Hn_V_S <- 1
    P_s_v <- x/sum(x)
    Hn_S_V <- entropy(P_s_v)
  } else {
    # Average mut.ind VOC
    { P_v = c()
    P_s_v = matrix(NA, NI, NJ)
    Hn_S_vj = c()
    P_v_s = matrix(NA, NI, NJ)
    Hn_V_si = c()
    for (j in 1: NJ) {
      P_v[j] <- colSums(x)[j]/sum(x)
      for (i in 1: NI) {
        P_s_v[i, j] <- x[i, j]/colSums(x)[j]
        P_v_s[i, j] <- x[i, j]/rowSums(x)[i]
      }
    }
    }
    P_s = 1/NI # assume equal abundance of each herbivore species
    Hn_S_vj <- apply(P_s_v, 2, entropy)  # scale matrix so that colsum is 1
    Hn_S_V <- sum(P_v*Hn_S_vj)
    Hn_V_si <- apply(P_v_s, 1, entropy)
    Hn_V_S  <- sum(P_s*Hn_V_si)
    Hn_V <- entropy(P_v)
  }
  results <- list(Hn_S_V, Hn_V_S, Hn_V)
  names(results) <- c("Hn_S_V", "Hn_V_S", "Hn_V")
  return(results)
}

# alternative function for conditional entropy

c_entropy <- function(mat,norm=T){
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

# calculate specific information

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


# calculate plant species fitness
IApV <- function(AP,AV, ii){
  A.p <- as.logical(AP[,ii]) #which herbivores are on this plant
  if (sum(A.p)==0){
    Fp <- 0
  } else{
    m<- nrow(AP)
    I <-specific_information(AV)
    Fp <-  1 - sum(sapply(1:m, function(x){
      (A.p[x]/sum(A.p))*I[[x]]})) # weighted average of specific information
  }
  return(Fp)
}

# figure label function
#https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

