#' Sci Data Paper Figure 5a/5b
#'
#' @import gridtext
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
    geom_hline(yintercept=-log10(0.05), col="green") + xlab(NULL) +
    ylab("-log(adjusted p-value)") + theme(legend.position="none",
                                       axis.text=element_text(size=rel(0.55)))
                                       #axis.text=element_text(size=10))

  # wrangle data for lollipop plot
  temp <- diffTab[diffTab$diffexpressed != "NO",]
  temp <- temp[order(temp$AutoDom.AD),]
  temp$Protein <- factor(temp$Protein, levels=temp$Protein)

  # lollipop plot
  p4 <- ggplot(temp, aes(x=Protein, y=AutoDom.AD,
                        label=round(AutoDom.AD, digits=2))) +
    geom_point(stat='identity', aes(col=diffexpressed), 
               #size=6) +
               size=2) +
    scale_color_manual(name="AutoDom ADD",
                       labels=c("Under Abundant", "Over Abundant"),
                       values=c("UP"="#DC3220", "DOWN"="#005AB5")) +
    geom_text(color="white", size=2) +
    labs(title= NULL, y=NULL, x="Protein Group",
         subtitle=NULL) +
    coord_flip() + theme_bw()+ theme(axis.text=element_text(size=rel(0.55)))
    #theme(axis.text.x = element_text(size = 10))

  ytext <- gridtext::richtext_grob("log2 Fold Change")

  p5 <- grid.arrange(p0, p4, ncol=2,   bottom=ytext)

  ggsave("fig/figure5a_b.pdf",p5, dpi=600, dev='pdf',
         height=9, width=15, units="in")
  dev.off()

  # warnings back on
  options(warn=w)
}
