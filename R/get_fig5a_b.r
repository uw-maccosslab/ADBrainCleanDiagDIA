#' Sci Data Paper Figure 5a/5b
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @return figure5a_b.pdf
#' @export
get_fig5a_b=function() {
  
  # warnings off
  w <- getOption("warn")
  options(warn=-1)
  
  if( !exists("fig"))
    dir.create("fig")
  
  # load data and filter significant protein groups
  data("diffTab", envir=environment())
  
  diffTab$Protein <- sapply(strsplit(rownames(diffTab), "\\|"), '[', 3)
  diffTab$Protein <- sapply(strsplit(diffTab$Protein, ",", fixed=T), '[', 1)
  diffTab$diffexpressed <- "NO"
  diffTab$diffexpressed[diffTab$AutoDom.AD > 1 &
                          diffTab$adj.P.Val < 0.05] <- "UP"
  diffTab$diffexpressed[diffTab$AutoDom.AD < -1 &
                          diffTab$adj.P.Val < 0.05] <- "DOWN"
  diffTab$Protein[diffTab$diffexpressed == "NO"] <- NA
  
  # volcano plot
  p0=ggplot(data=diffTab, aes(x=AutoDom.AD,
                              y=-log10(adj.P.Val),
                              col=diffexpressed,
                              label=Protein)) +
    geom_point() + theme_minimal() + ggrepel::geom_text_repel(size=3) +
    scale_color_manual(values=c("#005AB5", "black", "#DC3220")) +
    geom_vline(xintercept=c(-1, 1), col="green") +
    geom_hline(yintercept=-log10(0.05), col="green") + xlab("log2 Fold Change") +
    ylab("-log(adjusted p-value)") + theme(legend.position="none",
                                           axis.text=element_text(size=rel(0.55)))
  #axis.text=element_text(size=10))
  
  data("metaDT", envir=environment())
  metaDT <- metaDT[1:62,]
  rownames(metaDT) <- metaDT$Replicate
  idx <- grep("High|Sporadic", metaDT$Condition)
  metaDT <- metaDT[-idx,]
  data("noBatchProt", envir=environment())
  noBatchProt <- noBatchProt[,-idx]
  data("diffTab", envir=environment())
  
  # rename condition groups
  metaDT$Condition = gsub("Path", "ADNC", gsub("Control ", "HCF - ", 
                                               gsub("AD", "ADD", metaDT$Condition)))
  
  # filter data with fdr < 0.05 and logFC > 1
  temp <- diffTab[diffTab$adj.P.Val < 0.05,]
  temp <- temp[abs(temp$AutoDom.AD) > 1, ]
  sigNames <- rownames(temp)
  protSig <- noBatchProt[rownames(noBatchProt) %in% sigNames,
                         grep("TZR",colnames(noBatchProt))]
  rownames(protSig) <- sapply(strsplit(sapply(strsplit(rownames(protSig)," @ "),
                                              "[", 2), ","), "[", 1)
  rownames(protSig) <- sapply(strsplit(rownames(protSig), "|",fixed=T), "[", 3)
  
  # heatmap factors to be shown
  annotation <- data.frame(Condition=factor(metaDT$Condition),
                           `Cognitive Status`=factor(metaDT$`Cognitive Status`),
                           Age=metaDT$Age,
                           Sex=factor(metaDT$Sex),
                           Batch=factor(metaDT$Batch))
  rownames(annotation)=metaDT$Replicate
  heat_colors <- brewer.pal(11, "RdBu")
  annot_colors=list(Condition=c(`AutoDom ADD`="#ff4500",
                                `HCF - Low ADNC`="#1874cd"))
  
  p4 <- ggplotify::as.ggplot(pheatmap(protSig, color=heat_colors,
                             cluster_rows=T, show_rownames=T,
                             annotation_col=annotation,
                             annotation_colors=annot_colors,
                             #clustering_distance_rows="correlation",
                             border_color=NA, fontsize=7, scale="row",
                             fontsize_row=5, height=20, show_colnames=F))
  
  p5 = grid.arrange(p0,p4,ncol=2)
  

  
  ggsave("fig/figure5a_b.pdf",p5, dpi=600, dev='pdf',
         height=9, width=15, units="in")
  dev.off()
  
  # warnings back on
  options(warn=w)
}
