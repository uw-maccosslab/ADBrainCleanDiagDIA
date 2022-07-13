#' Sci Data Paper Figure 2
#' @import grid
#' @import gridExtra
#' @import cowplot
#' @import ggplot2
#' @import reshape2
#'
#' @return figure2.pdf
#'
#' @export
get_fig2=function() {

  # warnings off
  w <- getOption("warn")
  options(warn=-1)

  if( !exists("fig"))
    dir.create("fig")

  # load data into local environment
  data(raw, envir=environment())
  data(totalID, envir=environment())
  data(runOrder, envir=environment())
  data(raw, envir=environment())
  data(median, envir=environment())
  data(metaDT, envir=environment())

  # wrangle
  runOrder <- subset(runOrder, runOrder$`Sample Name` %in% colnames(raw))
  runOrder$`Run Order` <- 1:nrow(runOrder)
  colnames(runOrder)[colnames(runOrder) == "Sample Name"] <- "Replicate"

  metaDT <- metaDT[order(match(metaDT$Condition,
                              c("AutoDom AD",
                                "Sporadic AD",
                                "Control High Path",
                                "Control Low Path",
                                "SMTG Reference",
                                "Brain Reference"))),]
  metaDT$Batch <- paste0("Batch ", metaDT$Batch)

  raw <- as.matrix(raw)[,match(metaDT$Replicate, colnames(raw))]

  totalID$class <- "Run Level FDR"
  rownames(totalID) <- totalID$Replicate

  temp <- rbind(
    rbind(cbind(`Peptide Detections`=colSums(raw == 0),
                Replicate=colnames(raw), class="Zero"),
          cbind(`Peptide Detections`=colSums(raw != 0)-colSums(raw == 0),
                Replicate=colnames(raw), class="Non-zero"),
          cbind(`Peptide Detections`= colSums(raw == 0),
                Replicate=colnames(raw), class="Experiment Level FDR"))
  )

  temp <- merge(temp, metaDT[,c("Replicate", "Condition", "Batch")],
                by="Replicate")
  temp <- merge(temp, runOrder, by="Replicate")
  temp$`Peptide Detections` <- as.numeric(temp$`Peptide Detections`)
  temp$Replicate <- factor(temp$Replicate, levels=runOrder$Replicate)

  # === start p0 =====
  p0 <- ggplot2::ggplot(data=temp, aes(x=Replicate,
                                       y=`Peptide Detections`,
                                       fill=class)) +
    geom_bar(stat="identity")+facet_grid(~Batch, scales="free") +
    scale_fill_grey(start=0.8, end=0.2) + theme_bw() +
    theme(axis.text.x=element_text(size=rel(0.7),angle=45, hjust=1))

  totalID <- totalID[match(unique(temp$Replicate),totalID$Replicate),]
  errbar <- cbind(totalID$`Peptide Detections`, 'NA', 'NA')
  errbar <- as.numeric(t(errbar))

  p0 <- p0 + geom_errorbar(aes(y=errbar,
                              ymin=errbar,
                              ymax=errbar,  col="Run Level FDR"),
                          linetype=1, size=1)

  # === start p1 + p2 =====
  raw[raw < 1] <- NA
  raw <- as.matrix(log2(raw))[,sort(colnames(raw))]
  median <- as.matrix(median)[,sort(colnames(median))]
  metaDT <- metaDT[order(metaDT$Replicate),]

  metaDT$Condition <- gsub("Brain Reference", "Interexperiment QC",
                           metaDT$Condition)
  metaDT$Condition <- gsub("SMTG Reference", "Interbatch QC",
                           metaDT$Condition)

  new_colors <- cbind(Condition=c("AutoDom AD", "Sporadic AD",
                                  "Control High Path","Control Low Path",
                                  "Interexperiment QC", "Interbatch QC"),
                     Colors=c("#ff4500", "#ffa500",
                              "#000080", "#1874cd",
                              "#8c8c8c", "#cccccc"))

  raw_long <- reshape2::melt(raw, varnames=c("Peptide","Replicate"),
                             value.name="Abundance")
  raw_long <- merge(raw_long, metaDT, by="Replicate")
  raw_long$Replicate <- factor(raw_long$Replicate, levels=runOrder$Replicate)

  med_long <- reshape2::melt(median, varnames=c("Peptide","Replicate"),
                             value.name="Abundance")
  med_long <- merge(med_long, metaDT, by="Replicate")
  med_long$Replicate <- factor(med_long$Replicate, levels=runOrder$Replicate)

  #reorder legend
  raw_long$Condition <- factor(raw_long$Condition,
                              levels=unique(raw_long$Condition)[c(3,4,5,6,2,1)])
  med_long$Condition <- factor(med_long$Condition,
                              levels=unique(med_long$Condition)[c(3,4,5,6,2,1)])

  p1 <- ggplot(raw_long, aes(Replicate, Abundance, fill=Condition)) +
    geom_boxplot() + coord_flip() +
    cowplot::theme_cowplot(12)+ theme(axis.text.y=element_text(size=rel(0.75)))+
    theme(legend.position="none")+ ylab("Log2 Abundance") +
    scale_fill_manual(values=new_colors[,"Colors"])

  p2 <- ggplot(med_long, aes(Replicate, Abundance, fill=Condition)) +
    geom_boxplot() + coord_flip() +
    cowplot::theme_cowplot(12)+ theme(axis.text.y=element_text(size=rel(0.75)))+
    theme(legend.position="none")+ ylab("Log2 Abundance") +
    scale_fill_manual(values=new_colors[,"Colors"])

  p3 <- ggplot(med_long, aes(Replicate, Abundance, fill=Condition)) +
    geom_boxplot() + coord_flip() +
    cowplot::theme_cowplot(12)+ theme(axis.text.y=element_text(size=rel(0.75)))+
    theme(legend.position="left", legend.box="vertical") +
    guides(colour=guide_legend(nrow=1)) + ylab("Log2 Abundance") +
    scale_fill_manual(values=new_colors[,"Colors"])

  grobs <- ggplotGrob(p3)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

  # add legend
  bottom_row <- plot_grid(
    plot_grid(
      gridExtra::grid.arrange(
        ggplot(raw_long, aes(2^Abundance)) + geom_density() +
          cowplot::theme_cowplot(12)+ xlab("Abundance"),
        ggplot(raw_long, aes(Replicate, 2^Abundance)) + geom_boxplot() +
          coord_flip() + cowplot::theme_cowplot(12)+
          theme(axis.text.y=element_text(size=rel(0.75)))+
          theme(legend.position="none")+ ylab("Abundance")+
          scale_fill_manual(values=new_colors[,"Colors"]
                            [order(new_colors[,"Condition"])]),
        heights=c(1, 4),
        ncol =1,  top=grid::textGrob("Raw",gp=gpar(fontsize=12#,font=3
        ))),
      gridExtra::grid.arrange(
        ggplot(raw_long, aes(Abundance)) + geom_density() +
          cowplot::theme_cowplot(12) + xlab("Log2 Abundance"),
        p1,  heights=c(1, 4),
        ncol =1,  top=grid::textGrob("Log2 Transform",gp=gpar(fontsize=12#,font=3
        ))),
      gridExtra::grid.arrange(
        ggplot(med_long, aes(Abundance)) + geom_density() +
          cowplot::theme_cowplot(12)+ xlab("Log2 Abundance"),
        p2,  heights=c(1, 4),
        ncol=1,  top=grid::textGrob("Median Norm",gp=gpar(fontsize=12#,font=3
        ))),
      labels=c('B', '', ''), ncol=3, label_size=12,rel_widths=c(1,1,1)),
    legend, ncol=2, rel_widths=c(1, .1)
  )

  plot_grid(p0, bottom_row, labels=c('A' ,''),
            label_size=12, ncol=1, scale=0.975)

  ggsave("fig/figure2.pdf", dpi=600, dev='pdf', height=20, width=15, units="in")
  dev.off()

  # warnings back on
  options(warn=w)
}






