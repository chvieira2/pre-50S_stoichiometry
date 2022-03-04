setwd("K:/Collaborations/Rainer_Nikolay/20200211_Rainer_Mol_Cell/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(gplots)
library(data.table)
library(stringr)







#### proteinGroups table load and preparation ####

PG <- fread("../txt_MBR/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)

Ribo_ptns <- fread("../../List_Ribo_Ptns.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T, header = F)
Ribo_ptns <- as.data.frame(str_split_fixed(Ribo_ptns$V1, " ", 2))
names(Ribo_ptns) <- c("Gene.names", "Uniprot.ID")


PG$Gene.names <- sapply(strsplit(PG$Gene.names, ";"), "[", 1)
PG$Majority.protein.IDs <- sapply(strsplit(PG$Majority.protein.IDs, ";"), "[", 1)

PG <- PG[order(PG$Gene.names),]

IBAQ = c("iBAQ.1_1", "iBAQ.1_2", "iBAQ.1_3",
         "iBAQ.2_1", "iBAQ.2_2", "iBAQ.2_3",
         "iBAQ.3_1", "iBAQ.3_2", "iBAQ.3_3",
         "iBAQ.4_1", "iBAQ.4_2", "iBAQ.4_3",
         "iBAQ.5_1", "iBAQ.5_2", "iBAQ.5_3",
         "iBAQ.CTL_f8_1",
         "iBAQ.CTL_f8_2",
         "iBAQ.CTL_f8_3",
         "iBAQ.CTL_f9_1",
         "iBAQ.CTL_f9_2",
         "iBAQ.CTL_f9_3",
         "iBAQ.CTL_f10_1",
         "iBAQ.CTL_f10_2",
         "iBAQ.CTL_f10_3",
         "iBAQ.KORPLP_f6_1",
         "iBAQ.KORPLP_f6_2",
         "iBAQ.KORPLP_f6_3",
         "iBAQ.KORPLP_f7_1",
         "iBAQ.KORPLP_f7_2",
         "iBAQ.KORPLP_f7_3",
         "iBAQ.KORPLP_f8_1",
         "iBAQ.KORPLP_f8_2",
         "iBAQ.KORPLP_f8_3",
         "iBAQ.KORPLP_f9_1",
         "iBAQ.KORPLP_f9_2",
         "iBAQ.KORPLP_f9_3",
         "iBAQ.KORPLP_f10_1",
         "iBAQ.KORPLP_f10_2",
         "iBAQ.KORPLP_f10_3",
         "iBAQ.KORPLP_f11_1",
         "iBAQ.KORPLP_f11_2",
         "iBAQ.KORPLP_f11_3")



#### Filters ####
#filter out contaminants, Rev and only identified by site
PG <- subset(PG, Reverse != "+")
PG <- subset(PG, Potential.contaminant != "+")
PG <- subset(PG, Only.identified.by.site != "+")




#transform Intensity and ratios to Log2 (or 10 if you prefer)
PG[c(IBAQ)] = log2(PG[c(IBAQ)])
# change Inf and NaN values for NA
is.na(PG[c(IBAQ)]) <- sapply(PG[c(IBAQ)], is.infinite)

is.na(PG[c(IBAQ)]) <- sapply(PG[c(IBAQ)], is.nan)



# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(IBAQ)) {
  cat(IBAQ[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[IBAQ[i]]),]))
  cat("\t")
  cat(mean(PG[,IBAQ[i]], na.rm = T))
  cat("\n")
  rm(i)
}


#### Define PG_summ as working table ####

PG_summ <- subset(PG, select = c("Majority.protein.IDs",
                                 "Gene.names", IBAQ))






#### Normalization ####
# All FALSE leads to normalization to the median
{
  Ribo <- subset(PG_summ, Majority.protein.IDs %in% Ribo_ptns$Uniprot.ID)
  
  median_all <- median(as.matrix(Ribo[IBAQ]), na.rm = T)
  
  for (foo in IBAQ) {
    median_foo <- median(as.matrix(Ribo[foo]), na.rm = T)
    
    PG_summ[foo] <- PG_summ[foo] - (median_foo - median_all)
    
    rm(foo)
  }
  rm(Ribo, median_all)
}





#### Stochiometry calculation ####
PG_summ_rpl <- subset(PG_summ, Majority.protein.IDs %in% Ribo_ptns$Uniprot.ID)
stochiometry <- c()
for (foo in IBAQ) {
  median_foo <- median(as.matrix(PG_summ_rpl[foo]), na.rm = T)
  
  PG_summ[,paste0("Stochiometry.", foo)] <- log2(100) + (PG_summ[,foo] - median_foo)
  
  stochiometry <- c(stochiometry, paste0("Stochiometry.", foo))
  rm(foo, median_foo)
}

rm(PG_summ_rpl)


#### Calculations mean ####
stochiometry_mean <- c("Stochiometry.iBAQ.1",
                       "Stochiometry.iBAQ.2",
                       "Stochiometry.iBAQ.3",
                       "Stochiometry.iBAQ.4",
                       "Stochiometry.iBAQ.5",
                       "Stochiometry.iBAQ.CTL_f8", "Stochiometry.iBAQ.CTL_f9", "Stochiometry.iBAQ.CTL_f10",
                       "Stochiometry.iBAQ.KORPLP_f6", "Stochiometry.iBAQ.KORPLP_f7",
                       "Stochiometry.iBAQ.KORPLP_f8", "Stochiometry.iBAQ.KORPLP_f9",
                       "Stochiometry.iBAQ.KORPLP_f10", "Stochiometry.iBAQ.KORPLP_f11")
for (foo in stochiometry_mean) {
  PG_summ[foo] <- rowMeans(PG_summ[,stochiometry[grep(foo, stochiometry)]], na.rm = T)
  rm(foo)
}





stochiometry_mean_sd<- paste0("SD.", stochiometry_mean)

PG_summ[,stochiometry_mean_sd[1]] <- apply(PG_summ[,stochiometry[1:3]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[2]] <- apply(PG_summ[,stochiometry[4:6]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[3]] <- apply(PG_summ[,stochiometry[7:9]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[4]] <- apply(PG_summ[,stochiometry[10:12]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[5]] <- apply(PG_summ[,stochiometry[13:15]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[6]] <- apply(PG_summ[,stochiometry[16:18]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[7]] <- apply(PG_summ[,stochiometry[19:21]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[8]] <- apply(PG_summ[,stochiometry[22:24]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[9]] <- apply(PG_summ[,stochiometry[25:27]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[10]] <- apply(PG_summ[,stochiometry[28:30]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[11]] <- apply(PG_summ[,stochiometry[31:33]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[12]] <- apply(PG_summ[,stochiometry[34:36]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[13]] <- apply(PG_summ[,stochiometry[37:39]], 1, function(x) sd(x, na.rm = T))
PG_summ[,stochiometry_mean_sd[14]] <- apply(PG_summ[,stochiometry[40:42]], 1, function(x) sd(x, na.rm = T))




#### Defining proteins of interest ####
Proteins_of_interest_Uniprot.ID <- c(as.character(Ribo_ptns$Uniprot.ID),
                                     "P42641", #ObgE 
                                     "P0A8X0", #YjgA 
                                     "P33643", #RluD 
                                     "P0AAT6") #RsfS 

ptn_names <- c("rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV", "rplW", "rplX", "rplY", "rpmA", "rpmB", "rpmC", "rpmD", "rpmE", "rpmF", "rpmG", "rpmH", "rpmI", "rpmJ", "obgE", "rluD", "rsfS", "yjgA")



BiogenesisFactors <- fread("../../BiogenesisFactors.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)[,1]

Biogenesis_names <- PG_summ[PG_summ$Majority.protein.IDs %in% BiogenesisFactors, "Gene.names"]

PG_summ_Ribo <- subset(PG_summ, Majority.protein.IDs %in% c(Proteins_of_interest_Uniprot.ID,
                                                            BiogenesisFactors))



#### Plot functions ####
# Boxplots
MyBoxplot <- function(df, variables,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/15),
                      ptn = Proteins_of_interest_Uniprot.ID, ShowPtn = ifelse(!is.null(ptn), T, F),
                      Dots = T, Lines = T, Bars = F,
                      mainPlot = T, aspectRatio = 2.5,
                      InPercentage.y = F, sec_breaks.y = c(1,10,50, 100, 200), scaled_breaks.y = log2(sec_breaks.y)) {
  
  
  
  if (limits_breaks <= 0) limits_breaks = 1
  
  all <- melt(df,
              id.vars = c("Gene.names", "Majority.protein.IDs"), measure.vars=variables)
  
  all2 <- subset(all, Majority.protein.IDs %in% ptn)
  
  
  plot <- ggplot(all, aes(variable,value)) +
    
    ylab(AxisName_y) +
    coord_flip(ylim = c(limits_y_min, limits_y_max)) +
    scale_x_discrete(labels = LabelNames) +
    
    theme_bw() +
    theme(#panel.grid = element_blank(),
      plot.title = element_text(lineheight=.8,
                                #face="bold", 
                                vjust=0.5,
                                hjust = 0.5,
                                size=30),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(#face = "bold",
        color = "black",
        angle=0, 
        vjust=1,
        size=20),
      axis.title.x = element_text(#face="bold",
        size=25,
        hjust = 0.5,
        vjust = 1.5),
      axis.text.y  = element_text(#face = "bold",
        color = "black",
        angle=0, 
        vjust=0.4,
        size=20),
      aspect.ratio = aspectRatio)
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y) +
      theme(axis.text.x  = element_text(#face = "bold",
        color = "black",
        angle=45,
        hjust = .75,
        vjust=1,
        size=12))
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(limits_y_min, limits_y_max, limits_breaks))
    
  }
  
  
  if (mainPlot) {
    require(ggbeeswarm)
    plot <- plot +
      geom_quasirandom(dodge.width = 1, varwidth = TRUE, na.rm = T, size = 1,
                       alpha = .25, shape = 21, fill = rgb(0,0,0), color = rgb(0,0,0,0))
    #geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, binwidth = 0.5, na.rm = T) +
    
    #geom_jitter(fill = rgb(0,0,0,.25), na.rm = T, width = 0.3) +
    
    #geom_violin(fill = rgb(0,0,0,.25), na.rm = T) +
    
    #geom_boxplot(fill = rgb(0,0,0,0), width = 0.75, na.rm = T, notch = T) +
  }
  
  
  if (ShowPtn) {
    
    if (Dots) {
      plot <- plot +
        
        geom_point(data = all2,
                   aes(x = variable, col = Gene.names),
                   size = 1, na.rm = T)
    }
    
    if (Lines) {
      plot <- plot +
        geom_line(data = all2,
                  aes(x = variable, col = Gene.names, group = Gene.names),
                  size = 1, linetype = 1, na.rm = T, alpha = 0.5)
      
    }
    
  }
  
  
  
  
  return(plot)
  
}



# Barplots
MyBarplot <- function(df, variables,
                      AxisName_y = NULL, AxisName_x = NULL,
                      limits_y_min = ifelse(round(min(df[variables], na.rm = T)) <= 0, round(min(df[variables], na.rm = T)),0),
                      limits_y_max = ifelse(round(max(df[variables], na.rm = T)) >= 8, round(max(df[variables], na.rm = T)),8),
                      limits_breaks = round((limits_y_max+limits_y_min)/15),
                      ptn = Proteins_of_interest_Uniprot.ID, ShowPtn = ifelse(!is.null(ptn), T, F),
                      mainPlot = T, aspectRatio = 4,
                      ErrorBars = F, variables_SD = stochiometry_FC_sd,
                      InPercentage.y = F, sec_breaks.y = c(1,10,50, 100, 200), scaled_breaks.y = log2(sec_breaks.y)) {
  
  
  
  if (limits_breaks <= 0) limits_breaks = 1
  
  all <- melt(df,
              id.vars = c("Gene.names", "Majority.protein.IDs"),
              measure.vars = variables)
  
  all2 <- subset(all, Majority.protein.IDs %in% ptn)
  
  
  plot <- ggplot(all2,
                 aes(y=value, x=Gene.names)) +
    
    ylab(AxisName_y) +
    xlab(AxisName_x) +
    coord_flip(ylim = c(limits_y_min, limits_y_max)) +
    theme_bw() +
    theme(#panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(lineheight=.8,
                                #face="bold", 
                                vjust=0.5,
                                hjust = 0.5,
                                size=12),
      axis.title.y = element_text(size=12,
                                  hjust = 0.5,
                                  vjust = 0.5),
      axis.text.x  = element_text(#face = "bold",
        color = "black",
        angle=45,
        hjust = .75,
        vjust=1,
        size=12),
      axis.title.x = element_text(#face="bold",
        size=12,
        hjust = 0.5,
        vjust = 1.5),
      axis.text.y  = element_text(#face = "bold",
        color = "black",
        angle=0, 
        vjust=0.4,
        size=12),
      aspect.ratio = aspectRatio)
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y)
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(limits_y_min, limits_y_max, limits_breaks))
    
  }
  
  
  if (mainPlot) {
    plot <- plot + 
      
      geom_bar(data = all, position="dodge", stat="identity", aes(fill=variable))
  }
  
  
  if (ShowPtn) {
    plot <- plot +
      
      geom_bar(data = all2,
               position = "dodge",
               width = 0.75,
               stat = "identity", aes(fill=variable))
  }
  
  if (ErrorBars) {
    
    all_sd <- melt(df,
                   id.vars = c("Gene.names", "Majority.protein.IDs"),
                   measure.vars = variables_SD)
    names(all_sd) <- c("Gene.names","Majority.protein.IDs", "variable","SD")
    
    all_sd2 <- subset(all_sd, Majority.protein.IDs %in% ptn)
    
    
    all2 <- cbind(all2, all_sd2["SD"])
    
    all2$min <- all2$value - all2$SD
    all2$max <- all2$value + all2$SD
    
    plot <- plot +
      
      geom_errorbar(data = all2,
                    aes(ymin=min, ymax=max, group = variable),
                    width=0.5,
                    position=position_dodge(.75))
  }
  
  
  return(plot)
  
}



# Dotplots
MyDotplot <- function(df, variables,
                      AxisName_y = NULL, AxisName_x = NULL,
                      limits_y_min = ifelse(round(min(df[variables], na.rm = T)) <= 0, round(min(df[variables], na.rm = T)),0),
                      limits_y_max = ifelse(round(max(df[variables], na.rm = T)) >= 8, round(max(df[variables], na.rm = T)),8),
                      limits_breaks = round((limits_y_max+limits_y_min)/15),
                      ptn = ptn_names, ShowPtn = T, MeanLines = T,
                      aspectRatio = 0.25,
                      InPercentage.y = T, sec_breaks.y = c(1,2,4,8,16,32,64,128), scaled_breaks.y = log2(sec_breaks.y)) {
  
  
  
  if (limits_breaks <= 0) limits_breaks = 1
  
  all <- melt(df,
              id.vars = c("Gene.names", "Majority.protein.IDs"),
              measure.vars = variables)
  
  all2 <- subset(all, Gene.names %in% ptn)
  
  all2$Gene.names <- factor(all2$Gene.names, levels = ptn)
  all2$Sample <- factor(ifelse(grepl("CTL_f8", all2$variable), "CTL_f8",
                               ifelse(grepl("CTL_f9", all2$variable), "CTL_f9",
                                      ifelse(grepl("CTL_f10", all2$variable), "CTL_f10",
                                             ifelse(grepl("KORPLP_f6", all2$variable), "KO_f6",
                                                    ifelse(grepl("KORPLP_f7", all2$variable), "KO_f7",
                                                           ifelse(grepl("KORPLP_f8", all2$variable), "KO_f8",
                                                                  ifelse(grepl("KORPLP_f9", all2$variable), "KO_f9",
                                                                         ifelse(grepl("KORPLP_f10", all2$variable), "KO_f10",
                                                                                "KO_f11")))))))),
                        levels = c("CTL_f8", "CTL_f9", "CTL_f10",
                                   "KO_f6", "KO_f7", "KO_f8", "KO_f9", "KO_f10", "KO_f11"))
  
  all2$Replicate <- factor(ifelse(grepl("_1", all2$variable), "1",
                                  ifelse(grepl("_2", all2$variable), "2", "3")),
                           levels = c("1", "2", "3"))
  
 
  # Necessary for shifting dots position
  all2$Gene.names.ID <- apply(all2["Gene.names"], 1, function(x) grep(x, ptn))
  
  
  ###
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  
  all3 <- data_summary(all2, varname="value", 
                       groupnames=c("Gene.names", "Gene.names.ID", "Majority.protein.IDs", "Sample"))
  
  
  plot <- ggplot(all2,
                 aes(y=value, x=Gene.names.ID+0.5)) +
    geom_hline(yintercept = log2(100), linetype = "dashed") +
    
    ylab(AxisName_y) +
    xlab(AxisName_x) +
    coord_cartesian(ylim = c(limits_y_min, limits_y_max)) +
    scale_x_continuous(breaks = 1:length(ptn), labels = ptn,
                       expand = expansion(add = c(0, 0))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=12),
          axis.title.y = element_text(size=12,
                                      hjust = 0.5,
                                      vjust = 0.5),
          axis.text.x  = element_text(color = "black",
                                      angle=90,
                                      hjust = 1,
                                      vjust= 1,
                                      size=12),
          axis.title.x = element_text(size=12,
                                      hjust = 0.5,
                                      vjust = 2),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.4,
                                      size=12),
          aspect.ratio = aspectRatio)
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y)
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(limits_y_min, limits_y_max, limits_breaks))
    
  }
  
  
  if (ShowPtn) {
    plot <- plot +
      
      geom_point(aes(fill=Sample, group = Sample), shape = 21, color = "black",
                 position= position_dodge(width = .80, preserve = "total"),
                 na.rm = T, size = 2)
  }
  
  
  if (MeanLines) {
    plot <- plot +
      
      geom_crossbar(data = all3,
                 aes(y=value, x=Gene.names.ID+0.5, ymin = value, ymax = value,
                     #color=Sample,
                     group = Sample),
                 width = 2, na.rm = T,
                 position= position_dodge(width = .80, preserve = "total"))
  }
  
  return(plot)
  
}


#### Stochiometry dotplot ####
plot_ribo_CTL <- MyDotplot(df = PG_summ, variables = stochiometry[1:9], ptn = ptn_names,
                           AxisName_y = "Stochiometry (%)",
                           limits_y_min = -5,
                           limits_y_max = 10,
                           sec_breaks.y = c(1,5,10,20,50,100,200)) +
  theme(legend.position = "top")

plot_ribo_KO <- MyDotplot(df = PG_summ, variables = stochiometry[10:length(stochiometry)], ptn = ptn_names,
                           AxisName_y = "Stochiometry (%)",
                           limits_y_min = -5,
                           limits_y_max = 10,
                           sec_breaks.y = c(1,5,10,20,50,100,200)) +
  theme(legend.position = "top")

pdf(paste("Stochiometry_Ribo_Dotplot.pdf", sep = ""),
    useDingbats = F)
grid.arrange(plot_ribo_CTL, plot_ribo_KO,
             nrow = 2)
dev.off()




plot_biogenesis_CTL <- MyDotplot(df = PG_summ, variables = stochiometry[1:9],
                                 ptn = sort(unique(c(Biogenesis_names,"obgE", "rluD", "rsfS", "yjgA"))),
                                 AxisName_y = "Stochiometry (%)",
                                 limits_y_min = -13,
                                 limits_y_max = 10,
                                 sec_breaks.y = c(1,5,10,20,50,100,200)) +
  theme(legend.position = "top")

plot_biogenesis_KO <- MyDotplot(df = PG_summ, variables = stochiometry[10:length(stochiometry)],
                                 ptn = sort(unique(c(Biogenesis_names,"obgE", "rluD", "rsfS", "yjgA"))),
                                 AxisName_y = "Stochiometry (%)",
                                 limits_y_min = -13,
                                 limits_y_max = 10,
                                 sec_breaks.y = c(1,5,10,20,50,100,200)) +
  theme(legend.position = "top")



pdf(paste("Stochiometry_Biogenesis_Dotplot.pdf", sep = ""),
useDingbats = F)
grid.arrange(plot_biogenesis_CTL, plot_biogenesis_KO,
             nrow = 2)
dev.off()
rm(plot_ribo_CTL, plot_ribo_KO,
   plot_biogenesis_CTL, plot_biogenesis_KO)



#### Stochiometry heatmap ####
library("pheatmap")


# Ribosomal proteins
foo <- PG_summ[PG_summ$Gene.names %in% ptn_names, c("Gene.names", stochiometry_mean)]

foo<- foo[match(ptn_names, foo$Gene.names),]

foo[is.na(foo)] <- -Inf
mt <- 2**data.matrix(foo[stochiometry_mean])
row.names(mt) <- foo$Gene.names
mt <- t(mt)

pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         breaks = seq(0,200,2),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,200,50),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("RiboProteins_Heatmap.pdf", sep = ""))

pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         #breaks = seq(-2,2,.5),
         border_color = "white",
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         #legend = T, legend_breaks = seq(-2,2,.5),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("RiboProteins_Heatmap_scale.pdf", sep = ""))





# Biogenesis factors
all_names <- sort(unique(c(Biogenesis_names, "obgE", "rluD", "rsfS", "yjgA")))
foo <- PG_summ[PG_summ$Gene.names %in% all_names, c("Gene.names", stochiometry_mean)]

foo<- foo[match(all_names, foo$Gene.names),]

foo[is.na(foo)] <- -Inf
mt <- 2**data.matrix(foo[stochiometry_mean])
row.names(mt) <- foo$Gene.names
mt <- t(mt)

pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         breaks = seq(0,20,.2),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,20,4),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Biogenesis_Heatmap_20.pdf", sep = ""))


pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         breaks = seq(0,5,.05),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,5,1),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Biogenesis_Heatmap_5.pdf", sep = ""))


pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         breaks = seq(0,3,.03),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,3,.6),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Biogenesis_Heatmap_3.pdf", sep = ""))



pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 3,
         #breaks = seq(-2,2,.5),
         border_color = "white",
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         #legend = T, legend_breaks = seq(-2,2,.5),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Biogenesis_Heatmap_scale.pdf", sep = ""))


rm(mt, foo)


#### boxplot ####
plot <- MyBoxplot(df = PG_summ, variables = IBAQ, AxisName_y = "log2(Normalized IBAQ)",
                  LabelNames = IBAQ,
                  ShowPtn = F)


pdf(paste("Boxplot_IBAQ.pdf", sep = ""), useDingbats = F)
grid.arrange(plot,
             nrow = 1)
dev.off()
rm(plot)




PG_summ_filtered <- subset(PG_summ, grepl("rps", PG_summ$Gene.names))

plot <- MyBoxplot(df = PG_summ, variables = IBAQ, AxisName_y = "log2(Normalized IBAQ)",
                  LabelNames = IBAQ, ptn = PG_summ_filtered$Majority.protein.IDs,
                  ShowPtn = T , Dots = T, Lines = T)

pdf(paste("Boxplot_RPSs.pdf", sep = ""), useDingbats = F)
grid.arrange(plot,
             nrow = 1)
dev.off()
rm(plot)

PG_summ_filtered <- subset(PG_summ, (iBAQ.CTL_f8_1 >= 31) &
                             !(Majority.protein.IDs %in% Proteins_of_interest_Uniprot.ID) &
                             !grepl("rps", PG_summ$Gene.names))

plot <- MyBoxplot(df = PG_summ, variables = IBAQ, AxisName_y = "log2(Normalized IBAQ)",
                  LabelNames = IBAQ, ptn = PG_summ_filtered$Majority.protein.IDs,
                  ShowPtn = T , Dots = T, Lines = T)

pdf(paste("Boxplot_InterestingPTNs.pdf", sep = ""), useDingbats = F)
grid.arrange(plot,
             nrow = 1)
dev.off()
rm(plot)
rm(PG_summ_filtered)







#### Stochiometry barplot ####
plot <- MyBarplot(df = PG_summ_Ribo, variables = stochiometry_mean, ErrorBars = T, variables_SD = stochiometry_mean_sd,
                  AxisName_y = "Stoichiometry (%)", sec_breaks.y = c(1,5,10,20,50,100,200),
                  mainPlot = F, ShowPtn = T, InPercentage.y = T, aspectRatio = 2, limits_y_min = -1)

pdf(paste("Stochiometry_RiboProteins.pdf", sep = ""),
useDingbats = F)
grid.arrange(plot,
             nrow = 1)
dev.off()
rm(plot)




#### Individual Stochiometry plots ####
xfoo <- apply(PG_summ_Ribo["Majority.protein.IDs"], 1, function(x) {
  
  df = subset(PG_summ_Ribo, Majority.protein.IDs == x)
  
  MyBarplot(df = df, variables = stochiometry_mean, ErrorBars = T,
            variables_SD = stochiometry_mean_sd,
            AxisName_y = "Stochiometry (%)", sec_breaks.y = c(1,2,4,8,16,32,64,128),
            mainPlot = F, ShowPtn = T, InPercentage.y = T, aspectRatio = 0.75) +
    ggtitle(PG_summ_Ribo[PG_summ_Ribo$Majority.protein.IDs == x, "Gene.names"]) +
    theme(plot.title = element_text(hjust = 0.5))
  
})

pdf(paste("Stochiometry_RiboProteins_individual.pdf", sep = ""),
useDingbats = F, width = 5*2, height = 5*length(xfoo)/4)
do.call("grid.arrange", c(xfoo, ncol= 2))
dev.off()
rm(xfoo)



#### Write PG_summ ####
#PG_summ_foo <- subset(PG_summ, select = c("Majority.protein.IDs", "Gene.names", stochiometry_mean))
#PG_summ_foo[,stochiometry_mean] <- 2**PG_summ_foo[,stochiometry_mean]
fwrite(PG_summ, file = "PG_summ.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)
#rm(PG_summ_foo)
