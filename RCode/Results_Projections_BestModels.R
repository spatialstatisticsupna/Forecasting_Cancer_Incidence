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

## Number of areas
n<- length(unique(data$Geography_code))

model.type <- c("model1", "model2", "model3")
nt.type <- "t2"

library(ggplot2)
for (i in year.ahead) {
  pred.year <- (t.to+(i-2)):(t.last+(i-3))
  
  eval(parse(text = paste0("data.results.",model.type," <- c()")))
  IS <- c(NA,NA,NA)
  
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
    
    for (mt in 1:length(model.type)) {
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
      
      eval(parse(text = paste0("data2$Fitted.median.",model.type[mt]," <- res$shared_",nt.type,"$summary.fitted.values$`0.5quant`*10^5")))
      eval(parse(text = paste0("data2$Fitted.low.",model.type[mt]," <- res$shared_",nt.type,"$summary.fitted.values$`0.025quant`*10^5")))
      eval(parse(text = paste0("data2$Fitted.up.",model.type[mt]," <- res$shared_",nt.type,"$summary.fitted.values$`0.975quant`*10^5")))
      eval(parse(text = paste0("data2$Fitted.mean.",model.type[mt]," <- res$shared_",nt.type,"$summary.fitted.values$`mean`")))
      eval(parse(text = paste0("data2$Fitted.sd.",model.type[mt]," <- res$shared_",nt.type,"$summary.fitted.values$`sd`")))
      
      for (m in 1:n) {
        eval(parse(text = paste0("ARB <- abs((data2$Fitted.median.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]-data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)])/data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)])
        ")))
        eval(parse(text = paste0("IS_value <- data2$Fitted.up.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)] - data2$Fitted.low.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]")))
        # eval(parse(text = paste0("IS_",pred.year[l],"_",nt.type[ntt],"[1:length(period),k] <- rep(0,length(period))")))
        
        for (y in 1:i) {
          if(eval(parse(text = paste0("data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)]<data2$Fitted.low.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)]")))){
            eval(parse(text = paste0("IS[y]<- IS_value[y] + 2/0.05*(data2$Fitted.low.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)]-data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)])")))
          }
          else if(eval(parse(text = paste0("data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)]>data2$Fitted.up.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)]")))){
            eval(parse(text = paste0("IS[y]<- IS_value[y] + 2/0.05*(- data2$Fitted.up.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)] + data2$Observed.Rates[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-(i-y))] & data2$ID_type==1)])")))
          }
          else{
            IS[y]<- IS_value[y]
          }
          
        }
        
        eval(parse(text = paste0("mean.Fitted.Values <- data2$Fitted.mean.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]*data2$population[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]
        ")))
        eval(parse(text = paste0("sd.Fitted.Values <- sqrt(data2$population[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]*data2$Fitted.mean.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]+(data2$population[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]*data2$Fitted.sd.",model.type[mt],"[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)])^2)
        ")))
        eval(parse(text = paste0("DSS <- ((data2$counts[which(data2$Geography_code==n_code[m] & data2$Year%in%period[(length(period)-2):length(period)] & data2$ID_type==1)]-mean.Fitted.Values)/sd.Fitted.Values)^2+2*log(sd.Fitted.Values)")))
        
        
        eval(parse(text = paste0("data.results.",model.type[mt]," <- rbind(data.results.",model.type[mt],",
                                   cbind(rep(model.type[mt],i),
                                        period[(length(period)-2):length(period)],
                                        rep(pred.year[l],i),
                                        c(1:i),
                                        rep(n_code[m],i),
                                        ARB, IS, DSS))")))
      }
      
    }
    
    
  }
  
  data.results <- c()
  for (mt in 1:length(model.type)) {
    # eval(parse(text = paste0("data.results.",model.type[mt]," <- data.frame(data.results.",model.type[mt],")")))
    # eval(parse(text = paste0("colnames(data.results.",model.type[mt],")[1:4] <- c('model', 'Year', 'pred.year', 'Geography_code')")))
    eval(parse(text = paste0("data.results <- rbind(data.results,
                          data.results.",model.type[mt],")")))
  }
  data.results <- data.frame(data.results)
  colnames(data.results) <- c("model", "year", "pred.year", "year.ahead","Geography_code", "ARB", "IS", "DSS")
  
  
}


data.results$ARB <- as.numeric(data.results$ARB)
data.results$IS <- as.numeric(data.results$IS)
data.results$DSS <- as.numeric(data.results$DSS)

data <- aggregate(data.results[,c(6,7,8)] , FUN = "mean", by = list(data.results$model, 
                                                                    data.results$year.ahead,
                                                                    data.results$Geography_code))

break_points.MARB <- c(round(seq(round(min(data$ARB, na.rm = TRUE)-0.01,2), round(max(data$ARB, na.rm = TRUE)+0.01,2), length.out = 6),2))
# break_points.IS <- c(round(seq(round(min(df$IS)-0.01,2), round(max(df$IS)+0.01,2), length.out = 6),2))
break_points.IS <- round(as.vector(quantile(data$IS, probs = seq(0, 1, 1/5), na.rm = TRUE)),2)
break_points.IS[1] <- break_points.IS[1]-1
break_points.IS[length(break_points.IS)] <- break_points.IS[length(break_points.IS)]+1
break_points.DSS <- c(round(seq(round(min(data$DSS, na.rm = TRUE)-0.01,2), round(max(data$DSS, na.rm = TRUE)+0.01,2), length.out = 6),2))


################################################################################
###########                 Figure 4                ###########
################################################################################
## Projections for areas with missing values

load("../../Data/Information_missing_areas.Rdata")
data <- data.results
data$ARB[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- NA
data$IS[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- NA
data$DSS[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- NA



data$ARB <- as.numeric(data$ARB)
data$IS <- as.numeric(data$IS)
data$DSS <- as.numeric(data$DSS)

data.map <- aggregate(data[,c(6,7,8)] , FUN = "mean", by = list(data$model, 
                                                                data$year.ahead,
                                                                data$Geography_code))

colnames(data.map) <- c("model", "year.ahead","Geography_code", "ARB", "IS", "DSS")

library(dplyr)
df2<-left_join(data.map, carto, by = c('Geography_code' = 'Code'))

#df2 <- df2[,-c(1,2,3,4,7)]

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")


cols.MARB <- c("#ffffcc","#c2e699", "#78c679", "#31a354", "#006837")
cols.IS <- c("#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
cols.DSS <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
library(ggplot2)
year_hd.labs <- c("one-year-ahead", "two-years-ahead", "three-years-ahead")
names(year_hd.labs) <- c("1", "2", "3")
 

model.name <- unique(df$model)
for (mn in 1:length(model.name)) {
  s.p.3 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = ARB)) +
    theme_void()+
    labs(title = paste(model.name[mn]),
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 33, face = "bold",
                                hjust=0.5,color="gray35"),
      
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold", 
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 20, color="gray35"),
      legend.text = element_text(size = 15, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.MARB)),
                 guide = "bins",
                 breaks = break_points.MARB)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='ARB',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))
  
  
  eval(parse(text = paste0("ARB.",model.name[mn]," <- s.p.3")))
  
  
  s.p.2 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = IS))+
    theme_void()+
    labs(title = ' IS ',
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 45, face = "bold",
                                hjust=0.5,color="gray35"),
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=30,hjust=0.1,vjust=1, face = "bold",
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 25, color="gray35"),
      legend.text = element_text(size = 20, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.IS)),
                 guide = "bins",
                 breaks = break_points.IS)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='IS',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 75,
                                  barheight = 1))
  
  eval(parse(text = paste0("IS.",model.name[mn]," <- s.p.2")))
  
  s.p.1 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = DSS))+
    theme_void()+
    labs(title = 'DSS',
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 33, face = "bold",
                                hjust=0.5,color="gray35"),
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold",
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 20, color="gray35"),
      legend.text = element_text(size = 15, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.DSS)),
                 guide = "bins",
                 breaks = break_points.DSS)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='DSS',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))
  
  eval(parse(text = paste0("DSS.",model.name[mn]," <- s.p.1")))
  
}


library(ggpubr)
fig1 <- ggarrange(ARB.model1, ARB.model2, ARB.model3,
          nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")

ggplot2::ggsave(filename = paste0("./Figures/Figure4_ARB.pdf"),
                plot = fig1,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")



fig2 <- ggarrange(IS.model1, IS.model2, IS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure4_IS.pdf"),
                plot = fig2,
                device = "pdf",
                dpi = 600,
                width = 16,
                height = 24,
                units = "in")


fig3 <- ggarrange(DSS.model1, DSS.model2, DSS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure4_DSS.pdf"),
                plot = fig3,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")





################################################################################
###########                 Figure 5                ###########
################################################################################
## Projections for areas with  no missing values

load("../../Data/Information_missing_areas.Rdata")
data <- data.results
data$ARB[which(data$Geography_code%in%inf.unknown$unknown.code)] <- NA
data$IS[which(data$Geography_code%in%inf.unknown$unknown.code)] <- NA
data$DSS[which(data$Geography_code%in%inf.unknown$unknown.code)] <- NA



data$ARB <- as.numeric(data$ARB)
data$IS <- as.numeric(data$IS)
data$DSS <- as.numeric(data$DSS)

data.map <- aggregate(data[,c(6,7,8)] , FUN = "mean", by = list(data$model, 
                                                                data$year.ahead,
                                                                data$Geography_code))

colnames(data.map) <- c("model", "year.ahead","Geography_code", "ARB", "IS", "DSS")

library(dplyr)
df2<-left_join(data.map, carto, by = c('Geography_code' = 'Code'))

#df2 <- df2[,-c(1,2,3,4,7)]

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")

# break_points.MARB <- c(round(seq(round(min(df$ARB, na.rm = TRUE)-0.01,2), round(max(df$ARB, na.rm = TRUE)+0.01,2), length.out = 6),2))
# # break_points.IS <- c(round(seq(round(min(df$IS)-0.01,2), round(max(df$IS)+0.01,2), length.out = 6),2))
# break_points.IS <- round(as.vector(quantile(df$IS, probs = seq(0, 1, 1/5), na.rm = TRUE)),2)
# break_points.IS[1] <- break_points.IS[1]-1
# break_points.IS[length(break_points.IS)] <- break_points.IS[length(break_points.IS)]+1
# break_points.DSS <- c(round(seq(round(min(df$DSS, na.rm = TRUE)-0.01,2), round(max(df$DSS, na.rm = TRUE)+0.01,2), length.out = 6),2))


cols.MARB <- c("#ffffcc","#c2e699", "#78c679", "#31a354", "#006837")
cols.IS <- c("#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
cols.DSS <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
library(ggplot2)
year_hd.labs <- c("one-year-ahead", "two-years-ahead", "three-years-ahead")
names(year_hd.labs) <- c("1", "2", "3")



model.name <- unique(df$model)
for (mn in 1:length(model.name)) {
  s.p.3 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = ARB)) +
    theme_void()+
    labs(title = paste(model.name[mn]),
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 33, face = "bold",
                                hjust=0.5,color="gray35"),
      
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold", 
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 20, color="gray35"),
      legend.text = element_text(size = 15, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.MARB)),
                 guide = "bins",
                 breaks = break_points.MARB)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='ARB',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))
  
  
  eval(parse(text = paste0("ARB.",model.name[mn]," <- s.p.3")))
  
  
  s.p.2 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = IS))+
    theme_void()+
    labs(title = ' IS ',
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 45, face = "bold",
                                hjust=0.5,color="gray35"),
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=30,hjust=0.1,vjust=1, face = "bold",
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 25, color="gray35"),
      legend.text = element_text(size = 20, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.IS)),
                 guide = "bins",
                 breaks = break_points.IS)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='IS',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 75,
                                  barheight = 1))
  
  eval(parse(text = paste0("IS.",model.name[mn]," <- s.p.2")))
  
  s.p.1 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = DSS))+
    theme_void()+
    labs(title = 'DSS',
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 33, face = "bold",
                                hjust=0.5,color="gray35"),
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold",
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 20, color="gray35"),
      legend.text = element_text(size = 15, color="gray35"),
      legend.key = element_rect()) +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.DSS)),
                 guide = "bins",
                 breaks = break_points.DSS)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='DSS',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))
  
  eval(parse(text = paste0("DSS.",model.name[mn]," <- s.p.1")))
  
}


library(ggpubr)
fig1 <- ggarrange(ARB.model1, ARB.model2, ARB.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")

ggplot2::ggsave(filename = paste0("./Figures/Figure5_ARB.pdf"),
                plot = fig1,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")



fig2 <- ggarrange(IS.model1, IS.model2, IS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure5_IS.pdf"),
                plot = fig2,
                device = "pdf",
                dpi = 600,
                width = 16,
                height = 24,
                units = "in")


fig3 <- ggarrange(DSS.model1, DSS.model2, DSS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure5_DSS.pdf"),
                plot = fig3,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")


################################################################################
###########                 Figure Combined                ###########
################################################################################
## Projections for areas with missing values

load("../../Data/Information_missing_areas.Rdata")
data <- data.results
# data$ARB[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- data$ARB[which(!data$Geography_code%in%inf.unknown$unknown.code)]*10
# data$IS[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- NA
# data$DSS[which(!data$Geography_code%in%inf.unknown$unknown.code)] <- NA



data$ARB <- as.numeric(data$ARB)
# data$IS <- as.numeric(data$IS)
# data$DSS <- as.numeric(data$DSS)

data.map <- aggregate(data[,c(6,7,8)] , FUN = "mean", by = list(data$model, 
                                                                data$year.ahead,
                                                                data$Geography_code))


colnames(data.map) <- c("model", "year.ahead","Geography_code", "ARB", "IS", "DSS")

data.map$ARB[which(!data.map$Geography_code%in%inf.unknown$unknown.code)] <- data.map$ARB[which(!data.map$Geography_code%in%inf.unknown$unknown.code)]
data.map$ARB_M[which(data.map$Geography_code%in%inf.unknown$unknown.code)] <- data.map$ARB[which(data.map$Geography_code%in%inf.unknown$unknown.code)]

library(dplyr)
df2<-left_join(data.map, carto, by = c('Geography_code' = 'Code'))

#df2 <- df2[,-c(1,2,3,4,7)]

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")


cols.MARB <- c("#c2e699", "#78c679", "#31a354", "#006837","#003c30")
cols.MARB_M <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")


# cols.IS <- c("#c7eae5","#80cdc1", "#35978f", "#01665e", "#003c30")
# cols.DSS <- c("#fee08b","#fdae61", "#f46d43", "#d53e4f", "#9e0142")
library(ggplot2)
year_hd.labs <- c("one-year-ahead", "two-years-ahead", "three-years-ahead")
names(year_hd.labs) <- c("1", "2", "3")

library(ggnewscale)
model.name <- unique(df$model)
for (mn in 1:length(model.name)) {
  s.p.3 <- ggplot(subset(df, model == model.name[mn])) +
    geom_sf(aes(fill = ARB)) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.MARB)),
                 guide = "bins",
                 breaks = break_points.MARB)+
    new_scale_fill() +
    geom_sf(aes(fill = ARB_M)) +
    binned_scale(aesthetics = "fill", scale_name = "custom", 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.MARB_M)),
                 guide = "bins",
                 breaks = break_points.MARB)+
    # geom_polygon(aes_string(fill= "ARB_M"), size = 0.2)+
    # scale_fill_gradient(low ="pink", high ="red", na.value="blank") + 
    theme_void()+
    labs(title = paste(model.name[mn]),
         subtitle = "",
         x = NULL,
         y = NULL) +
    facet_wrap(~year_hd,
               labeller = labeller(year_hd =  year_hd.labs)) +
    # scale_fill_viridis(alpha=0.80)+
    theme(
      plot.title = element_text(size = 33, face = "bold",
                                hjust=0.5,color="gray35"),
      
      # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
      # plot.caption = element_text(size=10,color="gray35"),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold", 
                                  color="gray35"),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 20, color="gray35"),
      legend.text = element_text(size = 15, color="gray35"),
      legend.key = element_rect(),
      legend.box = "horizontal") +
    # scale_fill_binned( type = "viridis", breaks = break_points) +
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title='ARB',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title=' ',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=5,
                                  barwidth = 40,
                                  barheight = 1))
  
  
  eval(parse(text = paste0("ARB.",model.name[mn]," <- s.p.3")))
  
  
  # s.p.2 <- ggplot(subset(df, model == model.name[mn])) +
  #   geom_sf(aes(fill = IS))+
  #   theme_void()+
  #   labs(title = ' IS ',
  #        subtitle = "",
  #        x = NULL,
  #        y = NULL) +
  #   facet_wrap(~year_hd,
  #              labeller = labeller(year_hd =  year_hd.labs)) +
  #   # scale_fill_viridis(alpha=0.80)+
  #   theme(
  #     plot.title = element_text(size = 45, face = "bold",
  #                               hjust=0.5,color="gray35"),
  #     # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
  #     # plot.caption = element_text(size=10,color="gray35"),
  #     strip.text.x = element_text(size=30,hjust=0.1,vjust=1, face = "bold",
  #                                 color="gray35"),
  #     plot.background = element_rect(fill = "#f5f5f2", color = NA), 
  #     plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
  #     panel.background = element_rect(fill = "#f5f5f2", color = NA),
  #     panel.border = element_blank(),
  #     legend.position = 'bottom',
  #     legend.background = element_rect(fill = "#f5f5f2", color = NA),
  #     legend.title = element_text(size = 25, color="gray35"),
  #     legend.text = element_text(size = 20, color="gray35"),
  #     legend.key = element_rect()) +
  #   # scale_fill_binned( type = "viridis", breaks = break_points) +
  #   binned_scale(aesthetics = "fill", scale_name = "custom", 
  #                palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.IS)),
  #                guide = "bins",
  #                breaks = break_points.IS)+
  #   guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
  #                                 title='IS',  ##rename default legend
  #                                 title.position='top',
  #                                 title.hjust=0.5,
  #                                 ticks.colour='#f5f5f2',
  #                                 ticks.linewidth=5,
  #                                 barwidth = 75,
  #                                 barheight = 1))
  # 
  # eval(parse(text = paste0("IS.",model.name[mn]," <- s.p.2")))
  # 
  # s.p.1 <- ggplot(subset(df, model == model.name[mn])) +
  #   geom_sf(aes(fill = DSS))+
  #   theme_void()+
  #   labs(title = 'DSS',
  #        subtitle = "",
  #        x = NULL,
  #        y = NULL) +
  #   facet_wrap(~year_hd,
  #              labeller = labeller(year_hd =  year_hd.labs)) +
  #   # scale_fill_viridis(alpha=0.80)+
  #   theme(
  #     plot.title = element_text(size = 33, face = "bold",
  #                               hjust=0.5,color="gray35"),
  #     # plot.subtitle = element_text(size = 14,hjust=0.5,color="gray35"),
  #     # plot.caption = element_text(size=10,color="gray35"),
  #     strip.text.x = element_text(size=25,hjust=0.1,vjust=1, face = "bold",
  #                                 color="gray35"),
  #     plot.background = element_rect(fill = "#f5f5f2", color = NA), 
  #     plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
  #     panel.background = element_rect(fill = "#f5f5f2", color = NA),
  #     panel.border = element_blank(),
  #     legend.position = 'bottom',
  #     legend.background = element_rect(fill = "#f5f5f2", color = NA),
  #     legend.title = element_text(size = 20, color="gray35"),
  #     legend.text = element_text(size = 15, color="gray35"),
  #     legend.key = element_rect()) +
  #   # scale_fill_binned( type = "viridis", breaks = break_points) +
  #   binned_scale(aesthetics = "fill", scale_name = "custom", 
  #                palette = ggplot2:::pal_binned(scales::manual_pal(values = cols.DSS)),
  #                guide = "bins",
  #                breaks = break_points.DSS)+
  #   guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
  #                                 title='DSS',  ##rename default legend
  #                                 title.position='top',
  #                                 title.hjust=0.5,
  #                                 ticks.colour='#f5f5f2',
  #                                 ticks.linewidth=5,
  #                                 barwidth = 40,
  #                                 barheight = 1))
  # 
  # eval(parse(text = paste0("DSS.",model.name[mn]," <- s.p.1")))
  
}


library(ggpubr)
fig1 <- ggarrange(ARB.model1, ARB.model2, ARB.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")

ggplot2::ggsave(filename = paste0("./Figures/Figure_prueba2.pdf"),
                plot = fig1,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")



fig2 <- ggarrange(IS.model1, IS.model2, IS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure4_IS.pdf"),
                plot = fig2,
                device = "pdf",
                dpi = 600,
                width = 16,
                height = 24,
                units = "in")


fig3 <- ggarrange(DSS.model1, DSS.model2, DSS.model3,
                  nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
ggplot2::ggsave(filename = paste0("./Figures/Figure4_DSS.pdf"),
                plot = fig3,
                device = "pdf",
                dpi = 600,
                width = 15,
                height = 22,
                units = "in")




