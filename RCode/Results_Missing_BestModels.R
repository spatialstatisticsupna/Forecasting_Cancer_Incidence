
rm(list=ls())
library(ggplot2)
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

## Number of areas
n<- length(unique(data$Geography_code))


criterias <- c("MARB", "DSS", "IS")
model.type <- c("model1", "model2","model3")
nt.type <- "t2"

load("../../Data/Information_missing_areas.Rdata")


for (i in year.ahead) {
  pred.year <- (t.to+(i-2)):(t.last+(i-3))
  
  data.boxplot.MARB <- c()
  data.boxplot.DSS <- c()
  data.boxplot.IS <- c()
  for (mt in 1:length(model.type)) {
    eval(parse(text = paste0("MARBi_",nt.type," <- matrix(0, 
                               nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
    eval(parse(text = paste0("MRRMSEi_",nt.type," <- matrix(0, 
                               nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
    eval(parse(text = paste0("DSSi_",nt.type," <- matrix(0, 
                               nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
    eval(parse(text = paste0("ISi_",nt.type," <- matrix(0, 
                               nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
    eval(parse(text = paste0("CIi_",nt.type," <- matrix(0, 
                               nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
    
    
    
    for (l in 1:length(pred.year)) {
      eval(parse(text = paste0("MARB_",pred.year[l],"_",nt.type," <- matrix(0, 
                           nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
      eval(parse(text = paste0("MRRMSE_",pred.year[l],"_",nt.type," <- matrix(0, 
                           nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
      eval(parse(text = paste0("DSS_",pred.year[l],"_",nt.type," <- matrix(0, 
                           nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
      eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type," <- matrix(0, 
                           nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
      eval(parse(text = paste0("CI_",pred.year[l],"_",nt.type," <- matrix(0, 
                           nrow = length(t.from:(t.last-3)) , ncol = nrow(inf.unknown))")))
      
      
      period <- t.from:pred.year[l]
      data2 <- data[which(data$Year%in%period),]
      ## Number of areas
      n<- length(unique(data$Geography_code))
      
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
      data2$Observed.rates <- data2$counts/data2$population*10^5
      
      if(model.type[mt]=="model1"){
        eval(parse(text= paste0("load('./2_Multivariate_Models_indep/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_rw1.RData')")))
        res <- results.shared.unst.rw1
      }
      else if(model.type[mt]=="model2"){
        eval(parse(text= paste0("load('./2_Multivariate_Models_indep/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_scm.RData')")))
        res <- results.shared.unst.scm
      }
      else{
        eval(parse(text= paste0("load('./3_Multivariate_Models_shared_1/resul/",i,"_year_ahead/",t.from,"-",pred.year[l],"/results_shared_unst_rw1M.RData')")))
        res <- results.shared.unst.rw1M
      }
      
      for (ntt in 1:length(nt.type)) {
        eval(parse(text = paste0("data2$Fitted.median <- res$shared_",nt.type[ntt],"$summary.fitted.values$`0.5quant`*10^5")))
        eval(parse(text = paste0("data2$Fitted.mean <- res$shared_",nt.type[ntt],"$summary.fitted.values$mean")))
        eval(parse(text = paste0("data2$Fitted.sd <- res$shared_",nt.type[ntt],"$summary.fitted.values$sd")))
        eval(parse(text = paste0("data2$Fitted.low <- res$shared_",nt.type[ntt],"$summary.fitted.values$`0.025quant`*10^5")))
        eval(parse(text = paste0("data2$Fitted.up <- res$shared_",nt.type[ntt],"$summary.fitted.values$`0.975quant`*10^5")))
        
        period <- t.from:(pred.year[l]-i)
        data.missing <- data2[which(data2$Year%in%period & 
                                      data2$Geography_code%in%inf.unknown$unknown.code),] 
        
        
        for (k in 1:length(inf.unknown$unknown.code)) {
          mis.year <- 2001+(inf.unknown$V2[which(inf.unknown$unknown.code==inf.unknown$unknown.code[k])]-1)
          if(mis.year>(pred.year[l]-i)){
            mis.year <- (pred.year[l]-i)
          }
          
          eval(parse(text = paste0("MARB_",pred.year[l],"_",nt.type[ntt],"[1:min(inf.unknown$V2[k],length(period)), k] <- abs((data.missing$Fitted.median[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]-data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)])/data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)])")))
          eval(parse(text = paste0("MRRMSE_",pred.year[l],"_",nt.type[ntt],"[1:min(inf.unknown$V2[k],length(period)), k] <- ((data.missing$Fitted.median[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]-data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)])/data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)])^2")))
          
          IS_value <- data.missing$Fitted.up[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)] - data.missing$Fitted.low[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)]
          # eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type[ntt],"[1:length(period),k] <- rep(0,length(period))")))
          for (m in 1:min(inf.unknown$V2[k],length(period))) {
            if(data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]<data.missing$Fitted.low[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]){
              eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type[ntt],"[m,k]<- IS_value[m] + 2/0.05*(data.missing$Fitted.low[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)][m]-data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)][m])")))
            }
            else if(data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]>data.missing$Fitted.up[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]){
              eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type[ntt],"[m,k]<- IS_value[m] + 2/0.05*(- data.missing$Fitted.up[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)][m] + data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$ID_type==1)][m])")))
            }
            else{
              eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type[ntt],"[m,k]<- IS_value[m]")))
            }
            
            if(data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]>data.missing$Fitted.low[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)] & 
               data.missing$Observed.rates[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]<data.missing$Fitted.up[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year%in%period[m] & data.missing$ID_type==1)]){
              eval(parse(text = paste0("CI_",pred.year[l],"_",nt.type[ntt],"[m,k]<- 1")))
            }
          }
          
          mean.Fitted.Values <- data.missing$Fitted.mean[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]*data.missing$population[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]
          sd.Fitted.Values <- sqrt(data.missing$population[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]*data.missing$Fitted.mean[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]+(data.missing$population[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]*data.missing$Fitted.sd[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)])^2)
          eval(parse(text = paste0("DSS_",pred.year[l],"_",nt.type[ntt],"[1:min(inf.unknown$V2[k],length(period)), k] <- ((data.missing$counts[which(data.missing$Geography_code==inf.unknown$unknown.code[k] & data.missing$Year<=mis.year & data.missing$ID_type==1)]-mean.Fitted.Values)/sd.Fitted.Values)^2+2*log(sd.Fitted.Values)")))
        }
      }
    }
    
    groups.years <- rep(1:5, each =3)
    for (k in 1:length(inf.unknown$unknown.code)) {
      for (ntt in nt.type) {
        for (l in pred.year) {
          for (a in 1:(length(period)-1)) {
            if(eval(parse(text = paste0("MARB_",l,"_",ntt,"[a,k]")))!=0){
              data.boxplot.MARB <- rbind(data.boxplot.MARB,
                                         c(model.type[mt], k, l, 2000 + a,
                                           inf.unknown$V2[k],
                                           eval(parse(text = paste0("MARB_",l,"_",ntt,"[a,k]")))))
              
              data.boxplot.DSS <- rbind(data.boxplot.DSS,
                                        c(model.type[mt], k, l, 2000 + a,
                                          inf.unknown$V2[k],
                                          eval(parse(text = paste0("DSS_",l,"_",ntt,"[a,k]")))))
              
              
              data.boxplot.IS <- rbind(data.boxplot.IS,
                                       c(model.type[mt], k, l, 2000 + a,
                                         inf.unknown$V2[k],
                                         eval(parse(text = paste0("IS_",l,"_",ntt,"[a,k]")))))
            }
            
          }
        }
      }
    }
  }
  
  data.boxplot.MARB <- data.frame(data.boxplot.MARB)
  colnames(data.boxplot.MARB) <- c("model", "area", "validation.period", "year",
                                   "percentage.group", "value")
  data.boxplot.MARB$value <- as.numeric(data.boxplot.MARB$value)
  
  data.boxplot.DSS <- data.frame(data.boxplot.DSS)
  colnames(data.boxplot.DSS) <- c("model", "area", "validation.period", "year",
                                  "percentage.group", "value")
  data.boxplot.DSS$value <- as.numeric(data.boxplot.DSS$value)
  
  data.boxplot.IS <- data.frame(data.boxplot.IS)
  colnames(data.boxplot.IS) <- c("model", "area", "validation.period", "year",
                                 "percentage.group", "value")
  data.boxplot.IS$value <- as.numeric(data.boxplot.IS$value)
}


################################################################################
###########                 Figure 1                ###########
################################################################################
## Mean of the validation period for each area. One figure for each percentage of 
## missing data by area. In other words, one figure for each colouring of the areas
## We have 5 figures.

pd <- position_dodge(0.4)

# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)

for (cr in 1:length(criterias)) {
  eval(parse(text = paste0("data.",criterias[cr]," <- aggregate(data.boxplot.",criterias[cr],"$value,
  by = list(data.boxplot.",criterias[cr],"$model,
                                                            data.boxplot.",criterias[cr],"$area,
                                                            data.boxplot.",criterias[cr],"$year,
                                                            data.boxplot.",criterias[cr],"$percentage.group),
                         FUN = 'mean')")))
  
   
  eval(parse(text = paste0("colnames(data.",criterias[cr],") <- c('model', 'area',  'year', 'percentage.group', 'value')")))
  eval(parse(text = paste0("data.",criterias[cr],"$year <- factor(data.",criterias[cr],"$year)")))
 
  
  palet <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
             "#b3de69")
  # if (criteria[cr]=="MARB"){
  #   palet <- c("#d9ef8b","#a6d96a", "#66bd63", "#1a9850", "#006837", "#003c30", "#c7eae5")  
  # }
  # else if (criteria[cr]=="IS"){
  #   palet <- c("#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
  # }
  # else {
  #   palet <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
  # }
  eval(parse(text = paste0("p.group <- unique(data.",criterias[cr],"$percentage.group)")))
  
  eval(parse(text = paste0("max <- max(data.",criterias[cr],"$value)")))
  eval(parse(text = paste0("min <- min(data.",criterias[cr],"$value)")))
  eval(parse(text =paste0("pdf(file = 'Figures/Figure1_",criterias[cr],".pdf',
                          width = 25, height = 10)")))
  
  for (pg in p.group) {
    eval(parse(text = paste0("data <- subset(data.",criterias[cr],", percentage.group == ",pg,")")))
    areas <- unique(data$area)
    p1 <- ggplot() + ylim(min, max)
    for (a in 1:length(areas)) {
      p1 <- p1 + geom_point(data = subset(data, area == areas[a]), 
                            aes(x = year, y = value, colour = area, shape = model), 
                            position = pd,  size=7) 
    }
    
    p1 <- p1 + ggtitle(paste(pg,"years missing")) +
      labs(x = "Year", y = paste(criterias[cr]), colour = "Area ",
           shape = "Model Type") +
      scale_colour_manual(breaks = areas,
                          labels = areas,
                          values = palet) + 
      scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                         name=" ") + 
      # scale_shape_manual(labels = c(paste(selec.priors[o]," ",selec.method[o]), paste(model.SCM)), values = c(17,16)) +
      theme(text = element_text(size = 35), 
            legend.position = "top") 
    print(p1)
  }
  dev.off()
  
}


################################################################################
###########                 Figure 2                ###########
################################################################################
## Mean of the validation period for areas in each percentage of missing.  
##  In other words, mean for the areas for each colouring of the areas
## We have 1 figure.

pd <- position_dodge(0.4)


# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)

for (cr in 1:length(criterias)) {
  eval(parse(text = paste0("data.",criterias[cr]," <- aggregate(data.boxplot.",criterias[cr],"$value,
  by = list(data.boxplot.",criterias[cr],"$model,
                                                            
                                                            data.boxplot.",criterias[cr],"$year,
                                                            data.boxplot.",criterias[cr],"$percentage.group),
                         FUN = 'mean')")))
  
  
  eval(parse(text = paste0("colnames(data.",criterias[cr],") <- c('model',   'year', 'percentage.group', 'value')")))
  eval(parse(text = paste0("data.",criterias[cr],"$year <- factor(data.",criterias[cr],"$year)")))
  eval(parse(text = paste0("data.",criterias[cr],"$percentage.group <- factor(data.",criterias[cr],"$percentage.group)")))
  
  palet <- c("#990F0F", "#FF8F33", "#FFB2B2", "#A3CC51", "#51A3CC")
  # if (criterias[cr]=="MARB"){
  #   palet <- c("#d9ef8b","#a6d96a", "#66bd63", "#1a9850", "#006837")
  # }
  # else if (criterias[cr]=="IS"){
  #   palet <- c("#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
  # }
  # else {
  #   palet <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
  # }
  eval(parse(text = paste0("p.group <- unique(data.",criterias[cr],"$percentage.group)")))
  eval(parse(text = paste0("data <- data.",criterias[cr])))
  data <- data[order(data$percentage.group),]
  
  eval(parse(text =paste0("pdf(file = 'Figures/Figure2_",criterias[cr],".pdf',
                          width = 25, height = 10)")))
  
  p1 <- ggplot()
  for (pg in p.group) {
    p1 <- p1 + geom_point(data = subset(data, percentage.group == pg), 
                            aes(x = year, y = value, colour = percentage.group, shape = model), 
                            position = pd,  size=7) 
    }
    
    p1 <- p1 + ggtitle(" " ) +
      labs(x = "Year", y = paste(criterias[cr]), colour = "Years missing",
           shape = "Model Type") +
      scale_colour_manual(breaks = c("15", "12", "9", "6", "3"),
                          labels = c("15" = "15 years", "12" = "12 years", "9" = "9 years", 
                                     "6" = "6 years", "3" = "3 years"),
                          values = palet) + 
      scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                         name=" ") + 
      # scale_shape_manual(labels = c(paste(selec.priors[o]," ",selec.method[o]), paste(model.SCM)), values = c(17,16)) +
      theme(text = element_text(size = 35), legend.position = "top") 
    print(p1)
    dev.off()
  }


################################################################################
###########                 Figure 3                ###########
################################################################################
## One figure for each area.

pd <- position_dodge(0.4)


# pdf(file = paste0("./",cancer.location,"/figures/",removing.strategy[o],"_Assessment_Criteria.pdf"), width = 25 , height= 10)

for (cr in 1:length(criterias)) {
  eval(parse(text = paste0("data.",criterias[cr]," <- data.boxplot.",criterias[cr])))
  
  
  eval(parse(text = paste0("data.",criterias[cr],"$year <- factor(data.",criterias[cr],"$year)")))
  
  # palet <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  #            "#b3de69")
  if (criterias[cr]=="MARB"){
  palet <- c("#d9ef8b","#a6d96a", "#66bd63", "#1a9850", "#006837", "#003c30")
  }
  else if (criterias[cr]=="IS"){
    palet <- c("#f5f5f5","#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
  }
  else {
    palet <- c("#ffffbf","#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
  }
  
  eval(parse(text = paste0("areas <- unique(data.",criterias[cr],"$area)")))
  eval(parse(text =paste0("pdf(file = 'Figures/Figure3_",criterias[cr],".pdf',
                          width = 25, height = 10)")))
  
  for (a in 1:length(areas)) {
    eval(parse(text = paste0("data <- subset(data.",criterias[cr],", area == ",areas[a],")")))
    val.per <- unique(data$validation.period)
    
    p1 <- ggplot()
    for (vp in val.per) {
      p1 <- p1 + geom_point(data = subset(data, validation.period == vp), 
                            aes(x = year, y = value, colour = validation.period, shape = model), 
                            position = pd,  size=7) 
    }
    
    p1 <- p1 + ggtitle(paste("Area",areas[a])) +
      labs(x = "Year", y = paste(criterias[cr]), colour = "Validation Period") +
      scale_color_manual(labels = c("2001-2011", "2001-2012", "2001-2013", "2001-2014", 
                                    "2001-2015", "2001-2016"), values = palet)+
      scale_shape_manual(labels = c("Model 1","Model 2", "Model 3"), values=c(16, 17, 15),
                         name=" ") +
      # scale_shape_manual(labels = c(paste(selec.priors[o]," ",selec.method[o]), paste(model.SCM)), values = c(17,16)) +
      theme(text = element_text(size = 35), legend.position = "top") 
    print(p1)
  }
  dev.off()
  
}



