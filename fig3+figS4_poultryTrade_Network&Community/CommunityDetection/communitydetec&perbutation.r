library(igraph)

communitydetec_walktrap <- function(Gcsvfile,nstep=5){
  adjMatrixNet2<-as.matrix(read.csv(file=Gcsvfile,check.names=FALSE))
  n<-nrow(adjMatrixNet2)
  
  nb <-(n*(n-1)/2)
  m <- sum(adjMatrixNet2)/nb
  
  commdata <-list()
  commrawdata = NULL
  commnormal =NULL
  commdist_normal <- matrix(0,nrow = 100,ncol=5)
  # perbutation
  for (nnoise in 1:100){
    # generate Gaussian noise
    gnoise <- rnorm(nb, mean = 0, sd = 0.1*m)
    
    adjMnoised <- matrix(nrow=n,ncol=n)
    
    k=1
    
    for (j in 1:(n-1)){
      for (i in (j+1):n){
        # add Gaussian noise
        adjMnoised[i,j]<-adjMatrixNet2[i,j]+gnoise[k]
        k=k+1
      }
    }
    
    diag(adjMnoised) <- 0
    adjMnoised[adjMnoised<0]<-0
    # community detection from the network
    Net_normal <- graph_from_adjacency_matrix(adjMatrixNet2,mode="lower",weighted = TRUE)
    imc_normal <-cluster_walktrap(Net_normal,weights = E(Net_normal)$weight,steps = nstep)
    # community detection from perbutated networks
    Net2 <- graph_from_adjacency_matrix(adjMnoised, mode="lower",weighted = TRUE)
    imc2 <- cluster_walktrap(Net2,weights = E(Net2)$weight,steps = nstep)
    
    commnormal<-rbind(commnormal,c(nstep,as.vector(membership(imc_normal))))
    commrawdata <- rbind(commrawdata,c(nnoise,as.vector(membership(imc2))))
    # evaluation of community robustness
    comparemethods <- c("adjusted.rand","split.join","nmi","vi")
    
    for (comparemethod in comparemethods){
      commdist_normal[nnoise,comparemethod]<-compare(imc2,imc_normal,method = comparemethod)
      
    }
  }
  # results output
  write.table(file='commdist_normal_1213.csv',commdist_normal,row.names = TRUE,col.names = TRUE,sep=",")
  write.csv(commrawdata,'commrawdata_1213.csv',row.names = FALSE,col.names = FALSE,sep=",")
}
