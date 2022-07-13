#' Sci Data Paper Figure 7
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @return figure7.pdf
#' @export
get_fig7=function() {

  # warnings off
  w <- getOption("warn")
  options(warn=-1)

  if( !exists("fig"))
    dir.create("fig")

  # load data and wrangle
  data("metaDT", envir=environment())
  metaDT <- metaDT[1:62,]
  rownames(metaDT) <- metaDT$Replicate
  idx <- grep("High|Sporadic", metaDT$Condition)
  metaDT <- metaDT[-idx,]
  data("noBatchProt", envir=environment())
  noBatchProt <- noBatchProt[,-idx]
  data("diffTab", envir=environment())

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
  annot_colors=list(Condition=c(`AutoDom AD`="#ff4500",
                                `Control Low Path`="#1874cd"))

  hmap <- pheatmap(protSig, color=heat_colors,
           cluster_rows=T, show_rownames=T,
           annotation_col=annotation,
           annotation_colors=annot_colors,
           #clustering_distance_rows="correlation",
           border_color=NA, fontsize=7, scale="row",
           fontsize_row=5, height=20, show_colnames=F)


  ggsave("fig/figure7.pdf",hmap, dpi=600, dev='pdf',
         height=8, width=8, units="in")
  dev.off()

  # warnings back on
  options(warn=w)
}
