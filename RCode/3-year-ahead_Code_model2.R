################################################################################
###############       Spatio-temporal Shared models (INLA)       ###############
################################################################################


########################################
## graph for INLA         ##
########################################
library("spdep")
carto <- st_read("../Data/Carto_England/carto_england.shp")
carto <- carto[order(carto$Code), ]

## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
carto <- sf::st_as_sf(carto)

# Representamos y guardamos el grafo de vecindad en la carpeta de Figuras
W <- as(nb2mat(poly2nb(carto), style = "B"), "Matrix")

## Spatial neighborhood matrix (Q_{xi})
spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
g <- INLA::inla.read.graph("carto_nb.graph")
plot(g, "neato")

########################################
## Necessary constants         ##
########################################
t.from <- 2001 #first year
t.last <- 2019 #last year
t.to <- 2014 #first crossvalidation year

years.ahead <- 3

library("INLA")
inla.setOption(inla.mode = "classic")
repeat{
  ########################################
  ## Data organization for INLA         ##
  ########################################
  load("../../Data/Missing_Data/Lung_Male_Missing.Rdata")
  data <- data[order(data$Geography_code, data$Year),]
  period <- t.from:t.to
  data <- data[which(data$Year%in%period),]
  
  ## Number of areas
  n<- length(unique(data$Geography_code))
  
  ## Number of years
  t <- length(period)
  
  ##health-outcomes
  health <- c("Incidence_Count","Mortality_Count")
  J <- length(health)
  
  ##NA values for last year
  data[which(data$Year>(t.to-years.ahead)),"Incidence_Count"] <- NA
  
  counts <- c(data$Incidence_Count,data$Mortality_Count)
  health_outcome <- rep(health,each=n*t)
  data <- rbind(data,data)    
  data <- data[,-c(4,5)]
  data <- cbind(data,counts,health_outcome)
  
  ########################################
  ## Index for INLA         ##
  ########################################
  data$ID_type<-rep(c(1,2),each=n*t)   
  
  data$alpha1 <- rep(1:0,each=n*t)
  data$alpha2 <- rep(0:1,each=n*t)
  
  ID_area <- 1:n
  data$ID_area <- c(rep(ID_area,each = t),
                    rep(ID_area+n, each = t))
  data$ID_unst <- c(rep(NA,n*t),rep(ID_area,each = t))
  
  data$ID_year <- c(rep(1:t,n),rep(1:t,n))
  data$ID_year_scm <- c(rep(1:t,n),rep((t+1):(2*t),n))
  
  ## order
  data <- data[order(data[,'ID_type'], data[, 'ID_year'],data[, 'ID_area']),]
  
  ID_int <- as.vector(sapply(1:t, function(x) ID_area+(x-1)*n))
  data$ID_int1 <- c(ID_int,rep(NA,n*t))
  data$ID_int2 <- c(rep(NA,n*t),ID_int)
  
  
  
  ########################################
  ##                  ##
  ########################################
  ##Q
  Q_xi <- matrix(0, g$n, g$n)
  for (i in 1:g$n){
    Q_xi[i,i]=g$nnbs[[i]]
    Q_xi[i,g$nbs[[i]]]=-1
  }
  
  ##RW1
  D1 <- diff(diag(t), differences=1)
  Q_gammaRW1 <- t(D1)%*%D1
  
  ## constraints
  R_1_2 <- kronecker(Q_gammaRW1,diag(n))
  r_def_1_2 <- n
  A_constr_1_2<- kronecker(matrix(1,1,t), diag(n))
  A_constr_1_2 <- as(A_constr_1_2[-1,], "matrix")
  
  R_1_3 <- kronecker(diag(t),Q_xi)
  r_def_1_3 <- t
  A_constr_1_3<- kronecker(diag(t),matrix(1,1,n))
  A_constr_1_3 <- as(A_constr_1_3[-1,], "matrix")
  
  R_1_4 <- kronecker(Q_gammaRW1,Q_xi)
  r_def_1_4 <- n+t-1
  A.1.1 <- kronecker(matrix(1,1,t),diag(n))
  A.1.2 <- kronecker(diag(t),matrix(1,1,n))
  A_constr_1_4 <- as(rbind(A.1.1[-1,],A.1.2[-1,]), "matrix")
  
  
  #### Definimos los UNIFORM PRIORS
  sdunif="expression:

logdens=-log_precision/2;

return(logdens)"
  
  
  ###############################################################################
  #              MI6: SCM + unst mort + shared time + indep inter              ##
  ###############################################################################
  ########################################
  ##      Additive            ##
  ########################################
  formula = counts~ -1 + alpha1 + alpha2 +
    f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE) +
    f(ID_unst, model="iid", hyper=list(prec=list(prior=sdunif))) +
    f(ID_year_scm, model="besag2", graph=Q_gammaRW1, hyper=list(prec=list(prior=sdunif)),
    scale.model = T, constr = TRUE)

  shared_unst_scm_ad = inla(formula,
                            family = "poisson",
                            data = data,
                            E=population,
                            control.compute=list(dic=TRUE,
                                                 cpo=TRUE,
                                                 waic=TRUE,
                                                 hyperpar=TRUE,
                                                 config = TRUE,
                                                 return.marginals.predictor=TRUE),
                            control.inla=list(strategy="simplified.laplace",
                                              verbose = T,
                                              numint.maxfeval= 100000),
                            control.predictor = list(link=1,compute = TRUE))


  ########################################
  ##      Type I           ##
  ########################################
  formula = counts~ -1 + alpha1 + alpha2 +
    f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE) +
    f(ID_unst, model="iid", hyper=list(prec=list(prior=sdunif))) +
    f(ID_year_scm, model="besag2", graph=Q_gammaRW1, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE)  +
    f(ID_int1, model ="iid", constr = TRUE, hyper=list(prec=list(prior=sdunif))) +
    f(ID_int2, model ="iid", constr = TRUE, hyper=list(prec=list(prior=sdunif)))


  shared_unst_scm_t1 = inla(formula,
                            family = "poisson",
                            data = data,
                            E=population,
                            control.compute=list(dic=TRUE,
                                                 cpo=TRUE,
                                                 waic=TRUE,
                                                 hyperpar=TRUE,
                                                 config = TRUE,
                                                 return.marginals.predictor=TRUE),
                            control.inla=list(strategy="simplified.laplace",
                                              verbose = T,
                                              numint.maxfeval= 100000),
                            control.predictor = list(link=1,compute = TRUE))


  ########################################
  ##      Type II           ##
  ########################################
  formula = counts~ -1 + alpha1 + alpha2 +
    f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)), 
      scale.model = T, constr = TRUE) +
    f(ID_unst, model="iid", hyper=list(prec=list(prior=sdunif))) +
    f(ID_year_scm, model="besag2", graph=Q_gammaRW1, hyper=list(prec=list(prior=sdunif)), 
      scale.model = T, constr = TRUE)  +
    f(ID_int1, model="generic0", Cmatrix=R_1_2, rankdef=r_def_1_2, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_2, e=rep(0,n-1)))+
    f(ID_int2, model="generic0", Cmatrix=R_1_2, rankdef=r_def_1_2, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_2, e=rep(0,n-1)))
  
  
  shared_unst_scm_t2 = tryCatch({ inla(formula,
                                       family = "poisson",
                                       data = data,
                                       E=population,
                                       control.compute=list(dic=TRUE,
                                                            cpo=TRUE,
                                                            waic=TRUE,
                                                            hyperpar=TRUE,
                                                            config = TRUE,
                                                            return.marginals.predictor=TRUE),
                                       control.inla=list(strategy="simplified.laplace",
                                                         verbose = T,
                                                         numint.maxfeval= 100000),
                                       control.predictor = list(link=1,compute = TRUE))}, 
                                error = function(msg) {
                                  print(msg)
                                  return(NULL)
                                })
  
  
  
  ########################################
  ##      Type III           ##
  ########################################
  formula = counts~ -1 + alpha1 + alpha2 +
    f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE) +
    f(ID_unst, model="iid", hyper=list(prec=list(prior=sdunif))) +
    f(ID_year_scm, model="besag2", graph=Q_gammaRW1, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE)  +
    f(ID_int1, model="generic0", Cmatrix=R_1_3, rankdef=r_def_1_3, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_3, e=rep(0,t-1))) +
    f(ID_int2, model="generic0", Cmatrix=R_1_3, rankdef=r_def_1_3, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_3, e=rep(0,t-1)))


  shared_unst_scm_t3 = tryCatch({ inla(formula,
                                       family = "poisson",
                                       data = data,
                                       E=population,
                                       control.compute=list(dic=TRUE,
                                                            cpo=TRUE,
                                                            waic=TRUE,
                                                            hyperpar=TRUE,
                                                            config = TRUE,
                                                            return.marginals.predictor=TRUE),
                                       control.inla=list(strategy="simplified.laplace",
                                                         verbose = T,
                                                         numint.maxfeval= 100000),
                                       control.predictor = list(link=1,compute = TRUE))},
                                error = function(msg) {
                                  print(msg)
                                  return(NULL)
                                })


  ########################################
  ##      Type IV           ##
  ########################################
  formula = counts~ -1 + alpha1 + alpha2 +
    f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE) +
    f(ID_unst, model="iid", hyper=list(prec=list(prior=sdunif))) +
    f(ID_year_scm, model="besag2", graph=Q_gammaRW1, hyper=list(prec=list(prior=sdunif)),
      scale.model = T, constr = TRUE)  +
    f(ID_int1, model="generic0", Cmatrix=R_1_4, rankdef=r_def_1_4, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_4, e=rep(0,n+t-2))) +
    f(ID_int2, model="generic0", Cmatrix=R_1_4, rankdef=r_def_1_4, hyper=list(prec=list(prior=sdunif)),
      constr=TRUE, extraconstr=list(A=A_constr_1_4, e=rep(0,n+t-2)))


  shared_unst_scm_t4 = tryCatch({ inla(formula,
                                       family = "poisson",
                                       data = data,
                                       E=population,
                                       control.compute=list(dic=TRUE,
                                                            cpo=TRUE,
                                                            waic=TRUE,
                                                            hyperpar=TRUE,
                                                            config = TRUE,
                                                            return.marginals.predictor=TRUE),
                                       control.inla=list(strategy="simplified.laplace",
                                                         verbose = T,
                                                         numint.maxfeval= 100000),
                                       control.predictor = list(link=1,compute = TRUE))},
                                error = function(msg) {
                                  print(msg)
                                  return(NULL)
                                })



  
  ########################################
  ## save                               ##
  ########################################
  results.model2 <- list(shared_ad=shared_unst_scm_ad,
                                  shared_t1=shared_unst_scm_t1,
                                  shared_t2=shared_unst_scm_t2,
                                  shared_t3=shared_unst_scm_t3,
                                  shared_t4=shared_unst_scm_t4)
  save(results.model2, file=paste0("./resul/",years.ahead,"_year_ahead/",t.from,"-",t.to,"/results_model2.RData"))
  
  
  
  
  t.to <- t.to+1
  
  if(t.to==(t.last+1)){break}
}

