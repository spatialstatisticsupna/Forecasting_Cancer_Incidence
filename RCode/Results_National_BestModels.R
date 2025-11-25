rm(list=ls())

########################################
## Load Carto         ##
########################################
library("spdep")
carto <- st_read("../../Data/Carto_England/carto_england.shp")
carto <- carto[order(carto$Code), ]
carto <- sf::st_as_sf(carto)

########################################
## Necessary constants:         ##
########################################

t.from <- 2001 
t.to <- 2013
t.last <- 2019
K <- length(t.to:t.last)
K.years <- t.to:t.last

year.ahead <- 3#2:3

########################################
## Order data         ##
########################################
#Observed values
load("../../Data/Lung_Male_Aggregate.Rdata")
data <- data[order(data$Geography_code, data$Year),]
period <- t.from:t.last
data <- data[which(data$Year%in%period),]
nat.count <- aggregate(data$Incidence_Count, FUN = "sum", by = list(data$Year))

## Number of areas
n<- length(unique(data$Geography_code))

model.type <- c("model1","model2",
                "model3")
nt.type <- "t2"

for (i in year.ahead) {
  pred.year <- (t.to+(i-2)):(t.last+(i-3))
  ARB.data <- c()
  DSS.data <- c()
  CIL.data <- c()
  IS.data <- c()
  
  ARB.data.pro <- c()
  DSS.data.pro <- c()
  CIL.data.pro <- c()
  IS.data.pro <- c()
  for (l in 1:length(pred.year)) {
    
    period <- t.from:pred.year[l]
    data2 <- data[which(data$Year%in%period),]
    ## Number of areas
    n<- length(unique(data$Geography_code))
    n_code <- unique(data2$Geography_code)
    
    ## Number of years
    t <- length(period)
    
    ##health-outcomes
    health <- c("Incidence_Count","Mortality_Count")
    J <- length(health)
    
    
    counts <- c(data2$Incidence_Count,data2$Mortality_Count)
    health_outcome <- rep(health,each=n*t)
    data2 <- rbind(data2,data2)    
    data2 <- data2[,-c(4,5)]
    data2 <- cbind(data2,counts,health_outcome)
    
    data2$ID_type<-rep(c(1,2),each=n*t)
    data2 <- data2[order(data2$ID_type, data2$Year, data2$Geography_code),]
    data2$Observed.Rates <- data2$counts/data2$population*10^5
    
    national.table <- matrix(NA, nrow = 2*length(period), ncol = length(model.type)+2)
    cov.table <- matrix(0, nrow = length(period), ncol = length(model.type))
    
    for (mt in 1:length(model.type)) {
      if(model.type[mt]=="model1"){
        eval(parse(text= paste0("load('./2_Multivariate_Models_indep/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_rw1.RData')")))
        eval(parse(text = paste0("res <- results.shared.unst.rw1$shared_",nt.type)))
      }
      else if(model.type[mt]=="model2"){
        eval(parse(text= paste0("load('./2_Multivariate_Models_indep/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_scm.RData')")))
        eval(parse(text = paste0("res <- results.shared.unst.scm$shared_",nt.type)))
      }
      else{
        eval(parse(text= paste0("load('./3_Multivariate_Models_shared_1/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_rw1M.RData')")))
        eval(parse(text = paste0("res <- results.shared.unst.rw1M$shared_",nt.type)))
      }
      
      
      # start.time <- Sys.time()
      library(INLA)
      # set.seed(21032024)
      set.seed(12112024)
      n.samples <- 5000
      samples.rate <- inla.posterior.sample(n.samples, res)
      
      fun = function(...) {
        exp(Predictor)
      }
      samples.rate.exp = inla.posterior.sample.eval(fun, samples.rate)
      
      # set.seed(21032024)
      set.seed(12112024)
      eval(parse(text = paste0("count.",1:t," <- NULL")))
      for (a in 1:n) {
        for (y in 1:t) {
          eval(parse(text = paste0("count.",y,"[[",a,"]] <- rpois(n.samples,
                                   data2$population[which(data2$Geography_code == n_code[a] &
                                                      data2$Year == (2000+y) &
                                                      data2$ID_type == 1)]*samples.rate.exp[n*(y-1) + a,])")))
        }
      }
      
      ARB <- c()
      DSS <- c()
      CIL <- c()
      IS <- c()
      for (y in 1:t) {
        ##National incidence distribution:
        eval(parse(text = paste0("national.count <- Reduce(`+`,count.",y,")")))
        
        L <- quantile(national.count, probs = c(0.025,0.975))
        #HPD
        # library(coda)
        # L <- HPDinterval(as.mcmc(national.count), 0.95)
        
        national.table[2*(y-1)+1, 2*(mt-1) + 1] <- round(mean(national.count),0)
        national.table[2*(y-1)+2,2*(mt-1) + 1] <- paste("(",round(mean(L[1]),0)," - ",round(mean(L[2]),0),")")
        if(nat.count$x[y]>=mean(L[1]) & nat.count$x[y]<=mean(L[2])){
          cov.table[y,mt] <- 1
        }
        
        
        ARB[y] <- abs(round(mean(national.count),0)- nat.count$x[y])/nat.count$x[y]
        DSS[y] <- ((nat.count$x[y] - round(mean(national.count),0))/sd(national.count))^2 + 2*log(sd(national.count))
        CIL[y] <- round(mean(L[2]),0) - round(mean(L[1]),0)
        IS[y] <- round(mean(L[2]),0) - round(mean(L[1]),0)
        if(nat.count$x[y] < round(mean(L[1]),0)){
          IS[y] <- IS[y] + 2/0.025*(round(mean(L[1]),0) - nat.count$x[y])
        } else if(nat.count$x[y]> round(mean(L[2]),0)){
          IS[y] <- IS[y] + 2/0.025*(nat.count$x[y] - round(mean(L[2]),0))
        } 
      }
      
      # end.time <- Sys.time()
      # time.taken <- end.time - start.time
      # time.taken
      
      rm(list = c("samples.rate","samples.rate.exp","national.count"))
      
      ARB.data <- rbind(ARB.data, 
                        cbind(rep(model.type[mt], t), 2000+1:t, rep(l, t), ARB) )
      ARB.data.pro <- rbind(ARB.data.pro, ARB.data[which(ARB.data[,2]%in%c(period[t]-0:(i-1)) &
                                                           ARB.data[,3] == l),])
      ARB.data <- ARB.data[-which(ARB.data[,2]%in%c(period[t]-0:(i-1)) &
                                  ARB.data[,3] == l),]
  
      
      DSS.data <- rbind(DSS.data, 
                        cbind(rep(model.type[mt], t), 2000+1:t, rep(l, t), DSS) )
      DSS.data.pro <- rbind(DSS.data.pro,DSS.data[which(DSS.data[,2]%in%c(period[t]-0:(i-1)) &
                                                           DSS.data[,3] == l),] )
      DSS.data <- DSS.data[-which(DSS.data[,2]%in%c(period[t]-0:(i-1)) &
                                    DSS.data[,3] == l),]
      
      IS.data <- rbind(IS.data, 
                        cbind(rep(model.type[mt], t), 2000+1:t, rep(l, t), IS) )
      IS.data.pro <- rbind(IS.data.pro,IS.data[which(IS.data[,2]%in%c(period[t]-0:(i-1)) &
                                                       IS.data[,3] == l),] )
      IS.data <- IS.data[-which(IS.data[,2]%in%c(period[t]-0:(i-1)) &
                                    IS.data[,3] == l),]
      
      CIL.data <- rbind(CIL.data, 
                        cbind(rep(model.type[mt], t), 2000+1:t, rep(l, t), CIL) )
      CIL.data.pro <- rbind(CIL.data.pro, CIL.data[which(CIL.data[,2]%in%c(period[t]-0:(i-1)) &
                                                            CIL.data[,3] == l),])
      CIL.data <- CIL.data[-which(CIL.data[,2]%in%c(period[t]-0:(i-1)) &
                                    CIL.data[,3] == l),]
      
    }
    
    eval(parse(text = paste0("latex_table <- xtable::xtable(national.table, caption='Model selection criterias for  validation period ",l," ', digits=3)")))
    xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
  }

}


################################################################################
###########                 Table projections                ###########
################################################################################
ARB.data.pro <- data.frame(ARB.data.pro)
colnames(ARB.data.pro) <- c("model", "year", "validation.period", "value")
ARB.data.pro$value <- as.numeric(ARB.data.pro$value)

CIL.data.pro <- data.frame(CIL.data.pro)
colnames(CIL.data.pro) <- c("model", "year", "validation.period", "value")
CIL.data.pro$value <- as.numeric(CIL.data.pro$value)

DSS.data.pro <- data.frame(DSS.data.pro)
colnames(DSS.data.pro) <- c("model", "year", "validation.period", "value")
DSS.data.pro$value <- as.numeric(DSS.data.pro$value)

IS.data.pro <- data.frame(IS.data.pro)
colnames(IS.data.pro) <- c("model", "year", "validation.period", "value")
IS.data.pro$value <- as.numeric(IS.data.pro$value)


tab <- matrix(NA, ncol = 14, nrow = 3)
val.period <- 6
model.type <- c("model1", "model2", "model3")

for(mt in 1:length(model.type)){
  year1 <- 0
  year2 <- 0
  year3 <- 0
  
  CIL1 <- 0
  CIL2 <- 0
  CIL3 <- 0
  
  DSS1 <- 0
  DSS2 <- 0
  DSS3 <- 0
  
  IS1 <- 0
  IS2 <- 0
  IS3 <- 0
  
  res.model1 <- ARB.data.pro[which(ARB.data.pro$model == model.type[mt]),]
  CIL.model1 <- CIL.data.pro[which(CIL.data.pro$model == model.type[mt]),]
  DSS.model1 <- DSS.data.pro[which(DSS.data.pro$model == model.type[mt]),]
  IS.model1 <- IS.data.pro[which(IS.data.pro$model == model.type[mt]),]
  for (vp in 1:val.period) {
    year1 <- year1 + res.model1$value[which(res.model1$validation.period==vp)][1]
    year2 <- year2 + res.model1$value[which(res.model1$validation.period==vp)][2]
    year3 <- year3 + res.model1$value[which(res.model1$validation.period==vp)][3]
    
    CIL1 <- CIL1 + CIL.model1$value[which(CIL.model1$validation.period==vp)][1]
    CIL2 <- CIL2 + CIL.model1$value[which(CIL.model1$validation.period==vp)][2]
    CIL3 <- CIL3 + CIL.model1$value[which(CIL.model1$validation.period==vp)][3]
    
    DSS1 <- DSS1 + DSS.model1$value[which(DSS.model1$validation.period==vp)][1]
    DSS2 <- DSS2 + DSS.model1$value[which(DSS.model1$validation.period==vp)][2]
    DSS3 <- DSS3 + DSS.model1$value[which(DSS.model1$validation.period==vp)][3]
    
    IS1 <- IS1 + IS.model1$value[which(IS.model1$validation.period==vp)][1]
    IS2 <- IS2 + IS.model1$value[which(IS.model1$validation.period==vp)][2]
    IS3 <- IS3 + IS.model1$value[which(IS.model1$validation.period==vp)][3]
  }
  
  tab[1, 5*(mt-1) + 1] <- round(year1/val.period,3)
  tab[2, 5*(mt-1) + 1] <- round(year2/val.period,3)
  tab[3, 5*(mt-1) + 1] <- round(year3/val.period,3)
  
  tab[1, 5*(mt-1) + 2] <- CIL1/val.period
  tab[2, 5*(mt-1) + 2] <- CIL2/val.period
  tab[3, 5*(mt-1) + 2] <- CIL3/val.period
  
  tab[1, 5*(mt-1) + 3] <- round(DSS1/val.period,3)
  tab[2, 5*(mt-1) + 3] <- round(DSS2/val.period,3)
  tab[3, 5*(mt-1) + 3] <- round(DSS3/val.period,3)
  
  tab[1, 5*(mt-1) + 4] <- round(IS1/val.period,3)
  tab[2, 5*(mt-1) + 4] <- round(IS2/val.period,3)
  tab[3, 5*(mt-1) + 4] <- round(IS3/val.period,3)
  
}


eval(parse(text = paste0("latex_table <- xtable::xtable(tab, caption='Model selection criterias for  validation period ",l," ', digits=3)")))
xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))


################################################################################
###########                 Figure 6                ###########
################################################################################
library(ggplot2)
pd <- position_dodge(0.4)
data <- ARB.data
data <- data.frame(data)
colnames(data) <- c("model", "year", "validation.period", "value")
# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)
data$value <- as.numeric(data$value)

# palet <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
#            "#b3de69")
palet <- c("#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494")
 
p.group <- unique(data$validation.period)
  
max <- max(data$value)
min <- min(data$value)

pdf(file = 'Figures/Figure6.pdf', width = 25, height = 10)

p1 <- ggplot() + ylim(min, max) +
  geom_point(data = data, 
             aes(x = year, y = value, colour = validation.period, shape = model), 
             position = pd,  size=7) + 
  ggtitle("National Incidence") +
  labs(x = "Year", y = "ARB", colour = "Validation Period", shape = "Model Type") +
  scale_colour_manual(labels = c("2001-2011", "2001-2012", "2001-2013", "2001-2014", 
                                 "2001-2015", "2001-2016"), values = palet) + 
  scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                     name=" ") +
  theme(text = element_text(size = 35), 
        legend.position = "top")

print(p1)
 
dev.off()



################################################################################
###########                 Figure 7                ###########
################################################################################

pd <- position_dodge(0.4)
data <- DSS.data
data <- data.frame(data)
colnames(data) <- c("model", "year", "validation.period", "value")
# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)
data$value <- as.numeric(data$value)

# palet <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
#            "#b3de69")
palet <- c("#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494")

p.group <- unique(data$validation.period)

max <- max(data$value)
min <- min(data$value)

pdf(file = 'Figures/Figure7.pdf', width = 25, height = 10)

p1 <- ggplot() + ylim(min, max) +
  geom_point(data = data, 
             aes(x = year, y = value, colour = validation.period, shape = model), 
             position = pd,  size=7) + 
  ggtitle("National Incidence") +
  labs(x = "Year", y = "DSS", colour = "Validation Period", shape = "Model Type") +
  scale_colour_manual(labels = c("2001-2011", "2001-2012", "2001-2013", "2001-2014", 
                                 "2001-2015", "2001-2016"), values = palet) + 
  scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                     name=" ") +
  theme(text = element_text(size = 35), 
        legend.position = "top")

print(p1)

dev.off()



################################################################################
###########                 Figure 8                ###########
################################################################################

pd <- position_dodge(0.4)
data <- CIL.data
data <- data.frame(data)
colnames(data) <- c("model", "year", "validation.period", "value")
# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)
data$value <- as.numeric(data$value)

# palet <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
#            "#b3de69")
palet <- c("#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494")

p.group <- unique(data$validation.period)

max <- max(data$value)
min <- min(data$value)

pdf(file = 'Figures/Figure8.pdf', width = 25, height = 10)

p1 <- ggplot() + ylim(min, max) +
  geom_point(data = data, 
             aes(x = year, y = value, colour = validation.period, shape = model), 
             position = pd,  size=7) + 
  ggtitle("National Incidence") +
  labs(x = "Year", y = "CIL", colour = "Validation Period", shape = "Model Type") +
  scale_colour_manual(labels = c("2001-2011", "2001-2012", "2001-2013", "2001-2014", 
                                 "2001-2015", "2001-2016"), values = palet) + 
  scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                     name=" ") +
  theme(text = element_text(size = 35), 
        legend.position = "top")

print(p1)

dev.off()


