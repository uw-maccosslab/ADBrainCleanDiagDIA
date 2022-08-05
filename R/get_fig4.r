#' Sci Data Paper Figure 4
#'
#' @return figure4.pdf
#' @export
get_fig4 <- function() {

  # warnings off
  w <- getOption("warn")
  options(warn=-1)

  if( !exists("fig"))
    dir.create("fig")

  # load data
  data(raw, envir=environment())
  data(median, envir=environment())
  data(pepGrp, envir=environment())
  data(noBatchPep, envir=environment())
  data(protGrp, envir=environment())
  data(noBatchProt, envir=environment())
  data(metaDT, envir=environment())

  raw[log2(raw) <= 0] <- NA
  raw <- log2(raw)
  metaDT <- metaDT[order(metaDT$Replicate),]
  std.idx <- grep("RPR", metaDT$Replicate)

  #print("debug")

  raw <- as.matrix(raw)[,sort(colnames(raw))]
  median <- as.matrix(median)[,sort(colnames(median))]
  pepGrp <- as.matrix(pepGrp)[,sort(colnames(pepGrp))]
  noBatchPep <- as.matrix(noBatchPep)[,sort(colnames(noBatchPep))]
  protGrp <- as.matrix(protGrp)[,sort(colnames(protGrp))]
  noBatchProt <- as.matrix(noBatchProt)[,sort(colnames(noBatchProt))]

  raw[raw <= 0] <- NA
  pepGrp[pepGrp <= 0] <- NA
  noBatchPep[noBatchPep <= 0] <- NA
  protGrp[protGrp <= 0] <- NA
  noBatchProt[noBatchProt <= 0] <- NA

  # function to plot the line of best fit
  lowerFn <- function(data, mapping, method="lm", ...) {
    p <- ggplot2::ggplot(data, mapping) +
      ggplot2::geom_point(size=0.1, alpha=0.05) +
      ggplot2::geom_smooth(method=method, color="red", ...)
  }

  # scatter plot raw peptide control samples
  p1 <- GGally::ggpairs(data.frame(raw)[,std.idx],
               diag=list(continuous="blank"),
               upper=list(continuous=GGally::wrap("cor",fontface="bold",
                                                  size=5)),
               lower=list(continuous=GGally::wrap(lowerFn, method="lm"),
                            title="raw")) + ggplot2::theme_bw()

  # scatter plot batch adjusted peptide control samples
  p2=GGally::ggpairs(data.frame(noBatchPep)[,std.idx],
             diag=list(continuous="blank"),
             upper=list(continuous=GGally::wrap("cor",fontface="bold", size=5)),
             lower=list(continuous=GGally::wrap(lowerFn, method="lm"),
                          title="noBatch")) + ggplot2::theme_bw()

  g1 <- plot_grid(GGally::ggmatrix_gtable(p1),GGally::ggmatrix_gtable(p2),
                 ncol=2, labels=c('A' ,'B'), scale=0.975)

  ggplot2::ggsave("fig/figure4.pdf",g1, dpi=600, dev='pdf',
                  height=5, width=12, units="in")

  dev.off()

  # warnings back on
  options(warn=w)
}
