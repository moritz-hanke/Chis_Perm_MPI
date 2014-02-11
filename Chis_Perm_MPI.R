############ Parallelisungstests mit MPI

# Muss als Skript in der Form
# mpirun -n 1 R --slave -f Name_des_R-Skripts.R
# ausgeführt werden

# Der zu analysierende Datensatz sollte wie golgt aufgebaut sein:
#     - er heißt phaeno.snps
#     - es gibt eine Variable "status", die cases (1) & controls (0)
#       beinhaltet
#     - alle Variablen, die NICHT zu analysierende SNPs sind, belegen
#       Spalte 1 bis x; von Spalte x+1 bis length(phaeno.snps) [letzte
#       Spalte des Datensatzes] sind die SNPs

phaeno.snps <- read.csv("~/Arbeit/BIPS/GSEA/Test/SNPS.random", header=T, sep=";")

# first.snp gibt an, ab welcher Spalte die zu analysierende SNPs be-
# ginnen

first.snp <- 2

data.PER <- NULL   # Datensatz mit permutierten Sätzen
n.perm <- 5    #Anzahl Permutationen

namen <- names(phaeno.snps[,first.snp:(length(phaeno.snps)-900)])


library(parallel)

numWorkers <- 2



## Set up the 'cluster'
cl <- makeCluster(numWorkers, type = "MPI")

# ClusterExport muss gesetzt werden, damit die Variabeln in
# parLapply genutzt werden können
clusterExport(cl=cl, varlist=c("phaeno.snps", "first.snp",
                               "data.PER", "n.perm", "namen"), envir=environment())


# with a for-loop und parSapply
system.time(
chis.PER <-   parSapply(cl=cl, X=namen,
                        FUN=function(var){
                          
                          ps.perm <- vector(mode="integer", length=n.perm)
                          
                          
                          for(i in 1:n.perm){
                            phaeno.snps$status <- sample(phaeno.snps$status)
                            
                            if(length(levels(as.factor((phaeno.snps[,var])))) < 3){
                              #Chi^2 for Independence 2x3
                              chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                        simulate.p.value = T, B=2000)$p.value)
                              
                              smallest.p.per <- chi.2x3
                              names(smallest.p.per) <- var
                              smallest.p.per
                              
                            }else{
                              n10 <- sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="case"))
                              n00 <-sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="control"))
                              N.0 <- n10+n00
                              
                              n11 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="case"))
                              n01 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="control"))
                              N.1 <- n11+n01
                              
                              n12 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="case"))
                              n02 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="control"))
                              N.2 <- n12+n02
                              
                              N1 <- n10 + n11 + n12
                              N0 <- n00 + n01 + n02
                              N <- N1+ N0
                              
                              
                              # Statistken für rezessives, additives(codominant) und dominates
                              # Model von:
                              # http://scholarworks.gsu.edu/cgi/viewcontent.cgi?article=1073&context=math_theses
                              
                              # Trend-Test rezessiv
                              x <- c(0,0,1)
                              
                              rec <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                             N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                               N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                         (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.rec <- (1 - pchisq(rec, df=1))
                              
                              
                              # Trend-Test additiv/codominant
                              x <- c(0,1,2)
                              
                              codo <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                              N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                                N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                          (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.codo <- (1 - pchisq(codo, df=1))
                              
                              
                              # Trend-Test dominant
                              x <- c(0,1,1)
                              
                              dom <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                             N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                               N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                         (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.dom <- (1 - pchisq(dom, df=1))
                              
                              
                              #Chi^2 for Independence 2x3
                              
                              chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                        simulate.p.value = T, B=2000)$p.value)
                              
                              
                              
                              
                              smallest.p.per <- min(c(p.rec, p.codo, p.dom, chi.2x3))
                              names(smallest.p.per) <- var
                              smallest.p.per
                              
                            } #Ende else
                            
                            ps.perm[i] <- smallest.p.per
                            
                          } #Ende der Schleife
                          #names(ps.perm) <- var
                          ps.perm
                          
                        } # Ende der Funktion
                        
            , simplify = TRUE, USE.NAMES = TRUE
)  # Ende parLapply

)





# Observed Data mit parSapply
chis.OBS <-   parSapply(cl=cl, X=namen,
                        FUN=function(var){
                          
                            
                            if(length(levels(as.factor((phaeno.snps[,var])))) < 3){
                              #Chi^2 for Independence 2x3
                              chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                        simulate.p.value = T, B=2000)$p.value)
                              
                              smallest.p.per <- chi.2x3
                              names(smallest.p.per) <- var
                              smallest.p.per
                              
                            }else{
                              n10 <- sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="case"))
                              n00 <-sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="control"))
                              N.0 <- n10+n00
                              
                              n11 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="case"))
                              n01 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="control"))
                              N.1 <- n11+n01
                              
                              n12 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="case"))
                              n02 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="control"))
                              N.2 <- n12+n02
                              
                              N1 <- n10 + n11 + n12
                              N0 <- n00 + n01 + n02
                              N <- N1+ N0
                              
                              
                              # Statistken für rezessives, additives(codominant) und dominates
                              # Model von:
                              # http://scholarworks.gsu.edu/cgi/viewcontent.cgi?article=1073&context=math_theses
                              
                              # Trend-Test rezessiv
                              x <- c(0,0,1)
                              
                              rec <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                             N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                               N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                         (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.rec <- (1 - pchisq(rec, df=1))
                              
                              
                              # Trend-Test additiv/codominant
                              x <- c(0,1,2)
                              
                              codo <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                              N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                                N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                          (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.codo <- (1 - pchisq(codo, df=1))
                              
                              
                              # Trend-Test dominant
                              x <- c(0,1,1)
                              
                              dom <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                             N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                               N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                         (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                              
                              p.dom <- (1 - pchisq(dom, df=1))
                              
                              
                              #Chi^2 for Independence 2x3
                              
                              chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                        simulate.p.value = T, B=2000)$p.value)
                              
                              
                              
                              
                              smallest.p.per <- min(c(p.rec, p.codo, p.dom, chi.2x3))
                              smallest.p.per
                              
                            } #Ende else
                            
                          
                        } # Ende der Funktion
                        
                        , simplify = TRUE, USE.NAMES = TRUE
)  # Ende parLapply








# with a for-loop und parLapply
system.time(
  chis.PER <-   parLapply(cl=cl, X=namen,
                          fun=function(var){
                            
                            ps.perm <- vector(mode="integer", length=n.perm)
                            
                            
                            for(i in 1:n.perm){
                              phaeno.snps$status <- sample(phaeno.snps$status)
                              
                              if(length(levels(as.factor((phaeno.snps[,var])))) < 3){
                                #Chi^2 for Independence 2x3
                                chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                          simulate.p.value = T, B=2000)$p.value)
                                
                                smallest.p.per <- chi.2x3
                                names(smallest.p.per) <- var
                                smallest.p.per
                                
                              }else{
                                n10 <- sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="case"))
                                n00 <-sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="control"))
                                N.0 <- n10+n00
                                
                                n11 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="case"))
                                n01 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="control"))
                                N.1 <- n11+n01
                                
                                n12 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="case"))
                                n02 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="control"))
                                N.2 <- n12+n02
                                
                                N1 <- n10 + n11 + n12
                                N0 <- n00 + n01 + n02
                                N <- N1+ N0
                                
                                
                                # Statistken für rezessives, additives(codominant) und dominates
                                # Model von:
                                # http://scholarworks.gsu.edu/cgi/viewcontent.cgi?article=1073&context=math_theses
                                
                                # Trend-Test rezessiv
                                x <- c(0,0,1)
                                
                                rec <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                               N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                                 N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                           (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                                
                                p.rec <- (1 - pchisq(rec, df=1))
                                
                                
                                # Trend-Test additiv/codominant
                                x <- c(0,1,2)
                                
                                codo <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                                N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                                  N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                            (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                                
                                p.codo <- (1 - pchisq(codo, df=1))
                                
                                
                                # Trend-Test dominant
                                x <- c(0,1,1)
                                
                                dom <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                                               N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                                                 N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                                           (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
                                
                                p.dom <- (1 - pchisq(dom, df=1))
                                
                                
                                #Chi^2 for Independence 2x3
                                
                                chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                                                          simulate.p.value = T, B=2000)$p.value)
                                
                                
                                
                                
                                smallest.p.per <- min(c(p.rec, p.codo, p.dom, chi.2x3))
                                names(smallest.p.per) <- var
                                smallest.p.per
                                
                              } #Ende else
                              
                              ps.perm[i] <- smallest.p.per
                              
                            } #Ende der Schleife
                            #names(ps.perm) <- var
                            ps.perm
                            
                          } # Ende der Funktion
  )  # Ende parLapply
  
)




####################################################################
# replicate in parLapply


## Wo ist der erste SNP?
first.snp <- 2
## Namen der SNPs
namen <- names(phaeno.snps[,first.snp:length(phaeno.snps)])


## Anzahl an Workern 
numWorkers <- 2

## Set up the 'cluster'
cl <- makeCluster(numWorkers, type = "PSOCK")


################
#ä Hier der tricky part:
## Wrapping a function in a function:

## Funktion, die in parSapply mittels replicate x-mal ausgeführt werden soll:
chi_functions <-function(var){
  phaeno.snps$status <- sample(phaeno.snps$status)
  
  if(length(levels(as.factor((phaeno.snps[,var])))) < 3){
    #Chi^2 for Independence 2x3
    chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                              simulate.p.value = T, B=2000)$p.value)
    
    smallest.p.per <- chi.2x3
    names(smallest.p.per) <- var
    smallest.p.per
    
  }else{
    n10 <- sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="case"))
    n00 <-sum(na.omit(phaeno.snps[,var]==0 & phaeno.snps[,"status"]=="control"))
    N.0 <- n10+n00
    
    n11 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="case"))
    n01 <- sum(na.omit(phaeno.snps[,var]==1 & phaeno.snps[,"status"]=="control"))
    N.1 <- n11+n01
    
    n12 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="case"))
    n02 <- sum(na.omit(phaeno.snps[,var]==2 & phaeno.snps[,"status"]=="control"))
    N.2 <- n12+n02
    
    N1 <- n10 + n11 + n12
    N0 <- n00 + n01 + n02
    N <- N1+ N0
    
    
    # Statistken für rezessives, additives(codominant) und dominates
    # Model von:
    # http://scholarworks.gsu.edu/cgi/viewcontent.cgi?article=1073&context=math_theses
    
    # Trend-Test rezessiv
    x <- c(0,0,1)
    
    rec <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                   N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                     N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                               (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
    
    p.rec <- (1 - pchisq(rec, df=1))
    
    
    # Trend-Test additiv/codominant
    x <- c(0,1,2)
    
    codo <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                    N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                      N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                                (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
    
    p.codo <- (1 - pchisq(codo, df=1))
    
    
    # Trend-Test dominant
    x <- c(0,1,1)
    
    dom <- (N * (N*(n10*x[1]+n11*x[2]+n12*x[3]) -
                   N1*(N.0*x[1] + N.1*x[2] +N.2*x[3]))^2)/(
                     N1*N0*( N* (N.0*x[1]^2 + N.1*x[2]^2 +N.2*x[3]^2) -
                               (N.0*x[1]+N.1*x[2]+N.2*x[3])^2))
    
    p.dom <- (1 - pchisq(dom, df=1))
    
    
    #Chi^2 for Independence 2x3
    
    chi.2x3 <- try(chisq.test(xtabs(~ phaeno.snps$status + phaeno.snps[,var]),
                              simulate.p.value = T, B=2000)$p.value)
    
    
    
    
    smallest.p.per <- min(c(p.rec, p.codo, p.dom, chi.2x3))
    names(smallest.p.per) <- var
    smallest.p.per
    
  } #Ende else
}


# ClusterExport muss gesetzt werden, damit die Variabeln in
# parLapply genutzt werden können
clusterExport(cl=cl, varlist=c("phaeno.snps", "first.snp","namen",
                               "chi_functions"), envir=environment())


# with a for-loop und parSapply
system.time(
  chis.PER <-   parSapply(cl=cl, X=namen,
                          FUN=function(var){
                            replicate(50, chi_functions(var))
                          } # Ende der Funktion
                          
                          ,simplify = TRUE, USE.NAMES = TRUE
  )  # Ende parLapply
) 


## Shut down cluster
stopCluster(cl)


mpi.exit()
