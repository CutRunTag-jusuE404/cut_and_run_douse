# make figures

library(ggplot2)
library(viridis)


# SEACR_peaks -------------------------------------------------------------
rm(list=ls())
load("RData/02_SEACR_peaks.RData")

p1 = data.frame(ol = ctrl_olap.rate, samples = 1:length(ctrl_olap.rate))

plot1 =   ggplot(data=p1, aes(x=samples, y=ol/1000, group=1)) +
  geom_line(linetype = "dashed") +
  geom_point(color="red") +
  ylab("Number of peaks (x1,000)") +
  xlab("Samples") +
  labs(title = "Control samples D3 + WTT") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 1.5, y = p1$ol[1]/1000, label = as.character(p1$ol[1])) +
  annotate("text", x = 2.4, y = (p1$ol[2]/1000)+3, label = as.character(p1$ol[2])) +
  annotate("text", x = 3.3, y = (p1$ol[3]/1000)+4.5, label = as.character(p1$ol[3])) +
  annotate("text", x = 4.1, y = (p1$ol[4]/1000)+4, label = as.character(p1$ol[4])) +
  annotate("text", x = 5.1, y = (p1$ol[5]/1000)+4, label = as.character(p1$ol[5])) +
  annotate("text", x = 6, y = (p1$ol[6]/1000)+4, label = as.character(p1$ol[6])) +
  scale_x_discrete(name = "Samples", breaks = c(1:6), labels = c(1:6), limits = c(1:6))

tiff(filename = "figures/plot1_seacr_peaks_d3_wtt.tiff", units = "cm", width = 10, height = 10, res = 300)
plot1
dev.off()


p2 = data.frame(ol = d3_olap.rate, samples = 1:length(d3_olap.rate))

plot2 =   ggplot(data=p1, aes(x=samples, y=ol/1000, group=1)) +
  geom_line(linetype = "dashed") +
  geom_point(color="red") +
  ylab("Number of peaks (x1,000)") +
  xlab("Samples") +
  labs(title = "Control samples D3 peaks only") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 1.3, y = p1$ol[1]/1000, label = as.character(p1$ol[1])) +
  annotate("text", x = 2.2, y = (p1$ol[2]/1000)+3, label = as.character(p1$ol[2])) +
  annotate("text", x = 3.1, y = (p1$ol[3]/1000)+3, label = as.character(p1$ol[3])) +
  annotate("text", x = 3.95, y = (p1$ol[4]/1000)+3, label = as.character(p1$ol[4]))

tiff(filename = "figures/plot2_seacr_peaks_d3_only.tiff", units = "cm", width = 10, height = 10, res = 300)
plot2
dev.off()





# RUVseq ------------------------------------------------------------------

rm(list=ls())

library(RUVSeq)
library(cowplot)

load("RData/03_Ruvseq_diffEnrich.RData")


# normalisation -----------------------------------------------------------


set <- betweenLaneNormalization(cts.m, which = "upper")
s <- RUVSeq::RUVs(set, cIdx = rownames(set), scIdx = groups, k = k)

tiff(filename = "figures/plot3_ruvseq_NO_normalisation.tiff", units = "cm", width = 12, height = 12, res = 300)
plotPCA(set, col = colors[trt], main = "No Normalization PCA", labels = T)
dev.off()

tiff(filename = "figures/plot4_ruvseq_YES_normalisation.tiff", units = "cm", width = 12, height = 12, res = 300)
plotPCA(s$normalizedCounts, col = colors[trt], main = "Normalized PCA"
        , sub = paste("With", k ,"factors of unwanted variation removed", sep = " ")
        , labels = T, pch = 19)
dev.off()


# diff expression ---------------------------------------------------------


# plot5 -------------------------------------------------------------------


dt = res.tko.d3$table
dt$cutoff = cut(dt$logFC, c(-Inf, -1, 1, Inf), labels = c("< -1","-1 =< x =< 1","> 1"))
plot5 = ggplot(data = dt, mapping = aes(logCPM, logFC)) +
  geom_point(aes(colour = dt$cutoff), size = 0.6, alpha = 1/3) +
  scale_colour_discrete(name  = expression(paste("Log"[2],"FC cutoff"))) +
  labs(title = "Differential binding TKO/D3", x = "Relative Expression", y = expression(paste("Log"[2]," Fold Change"))) +
  theme(plot.title = element_text(hjust = 0.5))

tiff(filename = "figures/plot5_diff_binding_TKO_D3.tiff", units = "cm", width = 10, height = 10, res = 300)
plot5
dev.off()


# plot6 -------------------------------------------------------------------


dt = res.mut.d3$table
dt$cutoff = cut(dt$logFC, c(-Inf, -1, 1, Inf), labels = c("< -1","-1 =< x =< 1","> 1"))
plot6 = ggplot(data = dt, mapping = aes(logCPM, logFC)) +
  geom_point(aes(colour = dt$cutoff), size = 0.6, alpha = 1/3) +
  scale_colour_discrete(name  = expression(paste("Log"[2],"FC cutoff"))) +
  labs(title = "Differential binding mutT/D33", x = "Relative Expression", y = expression(paste("Log"[2]," Fold Change"))) +
  theme(plot.title = element_text(hjust = 0.5))

tiff(filename = "figures/plot6_diff_binding_mutT_D3.tiff", units = "cm", width = 10, height = 10, res = 300)
plot6
dev.off()



# plot7 -------------------------------------------------------------------

dt = res.mut.wtt$table
dt$cutoff = cut(dt$logFC, c(-Inf, -1, 1, Inf), labels = c("< -1","-1 =< x =< 1","> 1"))
plot7 = ggplot(data = dt, mapping = aes(logCPM, logFC)) +
  geom_point(aes(colour = dt$cutoff), size = 0.6) +
  scale_colour_discrete(name  = expression(paste("Log"[2],"FC cutoff"))) +
  labs(title = "Differential binding mutT/WTT", x = "Relative Expression", y = expression(paste("Log"[2]," Fold Change"))) +
  theme(plot.title = element_text(hjust = 0.5))

tiff(filename = "figures/plot7_diff_binding_mutT_WTT.tiff", units = "cm", width = 10, height = 10, res = 300)
plot7
dev.off()





# overlaps ----------------------------------------------------------------

load("RData/04_overlaps.RData")

# table1 ------------------------------------------------------------------

WriteXLS::WriteXLS(ExcelFileName = "tables/Table1_GenFeat_overlap.xls",OL1.table, col.names = T, row.names = F, AdjWidth = T)

# plot8 - barplot ---------------------------------------------------------

tiff(filename = "figures/plot8_barplot_GenFeat_overlap.tiff", units = "cm", width = 10, height = 10, res = 300)

ggplot(data = OL1.df, aes(x=factor(Group.1), y = x, fill = Group.1)) + 
  geom_bar(colour = "black", stat = "identity") +
  labs(y = "Number of genes", x = "Genomic feature") +
  theme(axis.title = element_text(size = 15, colour = "black")
        , axis.text = element_text(size = 11, colour = "black")
        , legend.position =  "none"
        #, legend.title = element_text(size = 15, face = "bold")
        #, legend.text = element_text(size = 13, face = "bold")
        #, legend.spacing.y = unit(0.1, "cm")
        #, legend.key.size = unit(0.8, "cm")
  ) +
  scale_fill_viridis(discrete = T, option = "B", direction = -1)

dev.off()





# table2  -----------------------------------------------------------------

WriteXLS::WriteXLS(ExcelFileName = "tables/Table2_protein_coding_genes_overlap.xls", pc.t , col.names = T, row.names = F, AdjWidth = T)


# table3 ------------------------------------------------------------------

WriteXLS::WriteXLS(ExcelFileName = "tables/Table3_lncRNA_genes_overlap.xls", linc.t , col.names = T, row.names = F, AdjWidth = T)




# plot9 - piechart1 -------------------------------------------------------
tiff(filename = "figures/plot9_piechart1_repeats.tiff", units = "cm", width = 10, height = 10, res = 300)
p1.plot.pie
dev.off()


# plot10 - piechart2 ------------------------------------------------------
tiff(filename = "figures/plot10_piechart2_repeats.tiff", units = "cm", width = 10, height = 10, res = 300)
p2.plot.pie
dev.off()


# gene expression ---------------------------------------------------------

load("RData/05_geneExpression.RData")

# plot11 ------------------------------------------------------------------
tiff(filename = "figures/plot11_scatterplot_RPKM.tiff", units = "cm", width = 15, height = 10, res = 300)

ggplot(fc1,aes(n,log.rpkm)) + geom_point(aes(color= "all genes" )) +
  geom_point(data=fc2,aes(color="genes on peaks")) +
  geom_hline(yintercept = median(fc1$log.rpkm), linetype="dashed", color = "blue", size = 1) +
  scale_color_manual(values = c("#999999", "#82206CFF")) +  # (viridis_pal()(2)) +
  labs(color="") +
  ylab(expression(paste("Log"[2],"RPKM"))) +
  xlab("Genes in alphabetical order") +
  theme(axis.text.x = element_blank()
        , axis.text.y = element_text(size = 10)
        , axis.ticks.x = element_blank()
        , axis.title = element_text(size = 10)
        , legend.text = element_text(size = 8)
  ) +
  annotate("text", label = "median RPKM", x = 10000, y = 3, color = "blue")
dev.off()