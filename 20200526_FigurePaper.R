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

PG_summ <- fread("PG_summ.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)

Ribo_ptns <- fread("../../List_Ribo_Ptns.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T, header = F)
Ribo_ptns <- as.data.frame(str_split_fixed(Ribo_ptns$V1, " ", 2))
names(Ribo_ptns) <- c("Gene.names", "Uniprot.ID")



IBAQ = c("iBAQ.1_1", "iBAQ.1_2", "iBAQ.1_3",
         "iBAQ.2_1", "iBAQ.2_2", "iBAQ.2_3",
         "iBAQ.3_1", "iBAQ.3_2", "iBAQ.3_3")#,
         #"iBAQ.4_1", "iBAQ.4_2", "iBAQ.4_3",
         #"iBAQ.5_1", "iBAQ.5_2", "iBAQ.5_3")
        # "iBAQ.CTL_f8_1",
        # "iBAQ.CTL_f8_2",
        # "iBAQ.CTL_f8_3",
        # "iBAQ.CTL_f9_1",
        # "iBAQ.CTL_f9_2",
        # "iBAQ.CTL_f9_3",
        # "iBAQ.CTL_f10_1",
        # "iBAQ.CTL_f10_2",
        # "iBAQ.CTL_f10_3",
        # "iBAQ.KORPLP_f6_1",
        # "iBAQ.KORPLP_f6_2",
        # "iBAQ.KORPLP_f6_3",
        # "iBAQ.KORPLP_f7_1",
        # "iBAQ.KORPLP_f7_2",
        # "iBAQ.KORPLP_f7_3",
        # "iBAQ.KORPLP_f8_1",
        # "iBAQ.KORPLP_f8_2",
        # "iBAQ.KORPLP_f8_3",
        # "iBAQ.KORPLP_f9_1",
        # "iBAQ.KORPLP_f9_2",
        # "iBAQ.KORPLP_f9_3",
        # "iBAQ.KORPLP_f10_1",
        # "iBAQ.KORPLP_f10_2",
        # "iBAQ.KORPLP_f10_3",
        # "iBAQ.KORPLP_f11_1",
        # "iBAQ.KORPLP_f11_2",
        # "iBAQ.KORPLP_f11_3")

stoichiometry <- c()
for (foo in IBAQ) {
  stoichiometry <- c(stoichiometry, paste0("Stochiometry.", foo))
  rm(foo)
}


PG_summ <- subset(PG_summ, select = c("Majority.protein.IDs", "Gene.names", 
                                      names(PG_summ)[names(PG_summ) %in% c(stoichiometry)]))







#### Stoichiometry ratio ####
ratio_stoichiometry <- c()
for (foo in IBAQ[1:6]) {
  ratio_stoichiometry <- c(ratio_stoichiometry, paste0("Ratio_Stoichiometry.", foo))
  

  
  rm(foo)
}


PG_summ[,ratio_stoichiometry[1]] <- (PG_summ[,stoichiometry[1]] - PG_summ[,stoichiometry[7]])+log2(100)
PG_summ[,ratio_stoichiometry[2]] <- (PG_summ[,stoichiometry[2]] - PG_summ[,stoichiometry[8]])+log2(100)
PG_summ[,ratio_stoichiometry[3]] <- (PG_summ[,stoichiometry[3]] - PG_summ[,stoichiometry[9]])+log2(100)
PG_summ[,ratio_stoichiometry[4]] <- (PG_summ[,stoichiometry[4]] - PG_summ[,stoichiometry[7]])+log2(100)
PG_summ[,ratio_stoichiometry[5]] <- (PG_summ[,stoichiometry[5]] - PG_summ[,stoichiometry[8]])+log2(100)
PG_summ[,ratio_stoichiometry[6]] <- (PG_summ[,stoichiometry[6]] - PG_summ[,stoichiometry[9]])+log2(100)









#### Calculations mean ####
stoichiometry_mean <- c("stoichiometry.iBAQ.1",
                        "stoichiometry.iBAQ.2",
                        "stoichiometry.iBAQ.3")#,
#"stoichiometry.iBAQ.4",
#"stoichiometry.iBAQ.5")

for (foo in stoichiometry_mean) {
  PG_summ[foo] <- rowMeans(PG_summ[,stoichiometry[grep(foo, stoichiometry)]], na.rm = T)
  rm(foo)
}




Ratio_stoichiometry_mean <- c("Ratio_Stoichiometry.iBAQ.1",
                              "Ratio_Stoichiometry.iBAQ.2")


PG_summ[,Ratio_stoichiometry_mean[1]] <- rowMeans(PG_summ[, ratio_stoichiometry[1:3]], na.rm = T)
PG_summ[,Ratio_stoichiometry_mean[2]] <- rowMeans(PG_summ[, ratio_stoichiometry[4:6]], na.rm = T)






#### Defining proteins of interest ####
Protein_interest_main <- c("P60438",#uL3(rpLE)
                           "P0ADY7",#uL16(rplP)
                           "P0A7M9",#bL31(rpmE)
                           "P0A7Q1",#bL35(rpmI)
                           "P0A7Q6",#bL36(rpmJ)
                           "P42641", #ObgE 
                           "P0AAT6", #RsfS
                           "P33643", #RluD
                           "P0A8X0", #YjgA
                           "P0AGK4") #YhbY)



Proteins_of_interest_Uniprot.ID <- c(as.character(Ribo_ptns$Uniprot.ID))#,
                                   #  "P42641", #ObgE 
                                    # "P0AAT6", #RsfS
                                    # "P33643", #RluD
                                   #  "P0A8X0", #YjgA
                                   #  "P0AGK4") #YhbY)

#ptn_names <- PG_summ[PG_summ$Majority.protein.IDs %in% Proteins_of_interest_Uniprot.ID, "Gene.names"]


ptn_names <- c("rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV", "rplW", "rplX", "rplY", "rpmA", "rpmB", "rpmC", "rpmD", "rpmE", "rpmF", "rpmG", "rpmH", "rpmI", "rpmJ")#,
             #  "obgE", "rsfS", "rluD", "yjgA", "yhbY")



BiogenesisFactors <- fread("../../BiogenesisFactors.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)[,1]

Biogenesis_names <- PG_summ[PG_summ$Majority.protein.IDs %in% BiogenesisFactors, "Gene.names"]

PG_summ_Ribo <- subset(PG_summ, Majority.protein.IDs %in% c(Proteins_of_interest_Uniprot.ID,
                                                            BiogenesisFactors))



#### Plot functions ####
# Dotplots
MyDotplot <- function(df, variables,
                      AxisName_y = "stoichiometry (%)", AxisName_x = NULL,
                      limits_y_min = ifelse(round(min(df[variables], na.rm = T)) <= 0, round(min(df[variables], na.rm = T)),0),
                      limits_y_max = ifelse(round(max(df[variables], na.rm = T)) >= 8, round(max(df[variables], na.rm = T)),8),
                      limits_breaks = round((limits_y_max+limits_y_min)/15),
                      ptn = ptn_names, ShowPtn = T, MeanLines = T, BackGroundBars = T,
                      aspectRatio = 0.25,
                      InPercentage.y = T, sec_breaks.y = c(1,2,4,8,16,32,64,128), scaled_breaks.y = log2(sec_breaks.y),
                      Dot.size = 3) {
  
  
  
  if (limits_breaks <= 0) limits_breaks = 1
  
  all <- melt(df,
              id.vars = c("Gene.names", "Majority.protein.IDs"),
              measure.vars = variables)
  
  all2 <- subset(all, Gene.names %in% ptn)
  
  all2$Gene.names <- factor(all2$Gene.names, levels = ptn)
  all2$Sample <- factor(ifelse(grepl(".iBAQ.1", all2$variable), "1",
                               ifelse(grepl(".iBAQ.2", all2$variable), "2",
                               ifelse(grepl(".iBAQ.3", all2$variable), "3",
                               ifelse(grepl(".iBAQ.4", all2$variable), "4",
                               ifelse(grepl(".iBAQ.5", all2$variable), "5",
                               ifelse(grepl("CTL_f8", all2$variable), "CTL_f8",
                               ifelse(grepl("CTL_f9", all2$variable), "CTL_f9",
                                      ifelse(grepl("CTL_f10", all2$variable), "CTL_f10",
                                             ifelse(grepl("KORPLP_f6", all2$variable), "KO_f6",
                                                    ifelse(grepl("KORPLP_f7", all2$variable), "KO_f7",
                                                           ifelse(grepl("KORPLP_f8", all2$variable), "KO_f8",
                                                                  ifelse(grepl("KORPLP_f9", all2$variable), "KO_f9",
                                                                         ifelse(grepl("KORPLP_f10", all2$variable), "KO_f10",
                                                                                "KO_f11"))))))))))))),
                        levels = c("1", "2", "3", "4", "5",
                                   "CTL_f8", "CTL_f9", "CTL_f10",
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
  
  
 { plot <- ggplot(all2,
                 aes(y=value, x=Gene.names.ID+0.5)) +
    geom_hline(yintercept = log2(100), linetype = "dashed") +
    
    ylab(AxisName_y) +
    xlab(AxisName_x) +
    coord_cartesian(ylim = c(limits_y_min, limits_y_max),
                    xlim = c(0.9, max(all2$Gene.names.ID, na.rm = T)+1.1)) +
    scale_x_continuous(breaks = 1:length(ptn), labels = ptn,
                       expand = expansion(add = c(0, 0))) +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
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
          aspect.ratio = aspectRatio)}
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y)
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(limits_y_min, limits_y_max, limits_breaks))
    
  }
  
  
  if (BackGroundBars) {
    
    
    
    for (i in seq(2,length(unique(all2$Gene.names)), 2)) {
      
      plot <- plot +
        
        annotate("rect", xmin=i, xmax=i+1, ymin=limits_y_min-1, ymax=limits_y_max+1, fill="black", alpha=0.2)
      
      rm(i)
    }
  }
  
  
  if (ShowPtn) {
    plot <- plot +
      
      geom_point(aes(fill=Sample, group = Sample), shape = 21, color = "black",
                 position= position_dodge(width = .80, preserve = "total"),
                 na.rm = T, size = Dot.size)
  }
  
  
  if (MeanLines) {
    plot <- plot +
      
      geom_crossbar(data = all3,
                 aes(y=value, x=Gene.names.ID+0.5, ymin = value, ymax = value,
                     #color=Sample,
                     group = Sample),
                 width = 1.5, na.rm = T, size = .5,
                 position= position_dodge(width = .80, preserve = "total"))
  }
  
  return(plot)
  
}


#### stoichiometry dotplot ####
plot_all_123 <- MyDotplot(df = PG_summ, variables = stoichiometry, ptn = ptn_names,
                          limits_y_min = -2,
                          limits_y_max = 9,
                          sec_breaks.y = c(0.1,1,5,20,50,100,200))

plot_all_12 <- MyDotplot(df = PG_summ, variables = stoichiometry[1:6], ptn = ptn_names,
                         limits_y_min = -2,
                         limits_y_max = 9,
                         sec_breaks.y = c(0.1,1,5,20,50,100,200))


pdf(paste("stoichiometry_RPLs.pdf", sep = ""),
    useDingbats = F)
grid.arrange(plot_all_123, plot_all_12,
             nrow = 2)
dev.off()






plot_ratio_12 <- MyDotplot(df = PG_summ, variables = ratio_stoichiometry, ptn = ptn_names,
                           limits_y_min = 3,
                           limits_y_max = 8,
                           sec_breaks.y = c(10,25,50,100,200),
                           AxisName_y = "stoichiometry ratio\nover mature 50S (%)")


pdf(paste("stoichiometry_ratio_RPLs.pdf", sep = ""),
    useDingbats = F)
grid.arrange(plot_all_123, plot_ratio_12,
             nrow = 2)
dev.off()
rm(plot_ratio_12, plot_all_12, plot_all_123)







plot_biogenesis_123 <- MyDotplot(df = PG_summ, variables = stoichiometry,
                                 ptn = sort(unique(c(Biogenesis_names,"obgE", "rluD", "rsfS", "yjgA"))),
                                 limits_y_min = -2,
                                 limits_y_max = 9,
                                 sec_breaks.y = c(0.1,1,5,20,50,100,200))

plot_biogenesis_12 <- MyDotplot(df = PG_summ, variables = stoichiometry[1:6],
                                 ptn = sort(unique(c(Biogenesis_names,"obgE", "rluD", "rsfS", "yjgA"))),
                                 limits_y_min = -2,
                                 limits_y_max = 9,
                                 sec_breaks.y = c(0.1,1,5,20,50,100,200))


pdf(paste("stoichiometry_Biogenesis.pdf", sep = ""),
useDingbats = F)
grid.arrange(plot_biogenesis_123, plot_biogenesis_12,
             nrow = 2)
dev.off()
rm(plot_biogenesis_123, plot_biogenesis_12,
   plot_biogenesis)



#### stoichiometry heatmap ####
library("pheatmap")


# 1 and 3
foo <- PG_summ[PG_summ$Majority.protein.IDs %in% Protein_interest_main, c("Gene.names", "Majority.protein.IDs", stoichiometry_mean[1:2])]

#Reorder
foo <- foo[match(Protein_interest_main, foo$Majority.protein.IDs),]

foo[is.na(foo)] <- -Inf

mt <- 2**data.matrix(foo[stoichiometry_mean[1:2]])
row.names(mt) <- foo$Gene.names
#mt <- t(mt)

pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 5,
         breaks = seq(0,120,1.2),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,120,30),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Paper_Heatmap_1x3.pdf", sep = ""))







# 5, 4 and 3
foo <- PG_summ[PG_summ$Majority.protein.IDs %in% Protein_interest_main, c("Gene.names", "Majority.protein.IDs", c("stoichiometry.iBAQ.5", "stoichiometry.iBAQ.4", "stoichiometry.iBAQ.3"))]

#Reorder
foo <- foo[match(Protein_interest_main, foo$Majority.protein.IDs),]

foo[is.na(foo)] <- -Inf

mt <- 2**data.matrix(foo[c("stoichiometry.iBAQ.5", "stoichiometry.iBAQ.4", "stoichiometry.iBAQ.3")])
row.names(mt) <- foo$Gene.names
#mt <- t(mt)

pheatmap(mt,
         #color = colorRampPalette(c("blue", "red"))(100),   #heatmap plotting colors
         #color = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)),
         color = rev(colorRampPalette(brewer.pal(10, "Spectral"))(100)),
         gaps_row = 5,
         breaks = seq(0,120,1.2),
         border_color = "white",
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         legend = T, legend_breaks = seq(0,120,30),
         cellwidth = 12.5, cellheight = 12.5,
         fontsize = 12,
         angle_col = "90",
         silent = F,
         filename = paste("Paper_Heatmap_5x4x3.pdf", sep = ""))