#' Sci Data Paper Figure 3
#'
#' @import plyr
#' @import grid
#' @import gridExtra
#' @import cowplot
#' @import ggplot2
#'
#' @return figure3.pdf
#' @export
get_fig3 <- function() {

  # warnings off
  w <- getOption("warn")
  options(warn=-1)

  if( !exists("fig") )
    dir.create("fig")

  # load data
  data(raw, envir=environment())
  data(median, envir=environment())
  data(pepGrp, envir=environment())
  data(noBatchPep, envir=environment())
  data(protGrp, envir=environment())
  data(noBatchProt, envir=environment())
  data(metaDT, envir=environment())

  metaDT <- metaDT[order(metaDT$Replicate),]
  std.idx <- grep("RPR", metaDT$Replicate)

  raw[raw < 1] <- NA
  pepGrp[pepGrp < 1] <- NA
  noBatchPep[noBatchPep < 1] <- NA
  protGrp[protGrp < 1] <- NA
  noBatchProt[noBatchProt < 1] <- NA

  raw <- as.matrix(log2(raw))[,sort(colnames(raw))]
  median <- as.matrix(median)[,sort(colnames(median))]
  pepGrp <- as.matrix(pepGrp)[,sort(colnames(pepGrp))]
  noBatchPep <- as.matrix(noBatchPep)[,sort(colnames(noBatchPep))]
  protGrp <- as.matrix(protGrp)[,sort(colnames(protGrp))]
  noBatchProt <- as.matrix(noBatchProt)[,sort(colnames(noBatchProt))]

  # format cv table
  get_cvTab <- function(mat, idx, text) {
    return(
      as.data.frame(cbind(`Log2 Median`=as.numeric(apply(
        mat[,idx], 1, median, na.rm=T)),
        CV=get_cv(as.matrix(2^mat[,idx])),
        method=rep(text,nrow(mat))))
    )
  }

  cv1 <- get_cvTab(raw, std.idx, "Raw")
  cv2 <- get_cvTab(median, std.idx, "Median Normalized")
  cv3 <- get_cvTab(pepGrp, std.idx, "PepGroup")
  cv4 <- get_cvTab(noBatchPep, std.idx, "Batch Adjusted")
  cv5 <- get_cvTab(protGrp, std.idx, "Protein Group")
  cv6 <- get_cvTab(noBatchProt, std.idx, "Batch Adjusted")

  df <- rbind(cv1,cv2,cv4)
  df$method <- factor(df$method, levels=c("Raw",
                                           "Median Normalized",
                                           "Batch Adjusted"))

  df[,1] <- as.numeric(df[,1])
  df[,2] <- as.numeric(df[,2])

  df <- df[-which(is.na(df[,1])),]
  df <- df[-which(is.na(df[,2])),]

  # common theme
  commonTheme <- list(labs(color="Density",fill="Density",
                          x="Log2 Median",
                          y="Coefficient of Variation",
                          title=NULL
  ),
  theme_bw(),
  theme(legend.position=c(0.995,0.995),
        legend.justification=c(0.995,0.995)))

  # peptide-centric density plot
  p1=ggplot(data=df,aes(`Log2 Median`,CV)) +
    geom_point(shape=16, size=0.1, color="grey66", show.legend=FALSE) +
    stat_density2d(aes(fill=..level..,alpha=..level..),
                   geom='polygon',
                   colour='black') +
    scale_fill_continuous(low="green",high="red")+
    geom_smooth(method="gam", formula=y ~ s(x, k=2),
                linetype=2,colour="red",se=F) +
    guides(alpha="none") +
    facet_grid(method~., scales="free_x") +
    commonTheme +
    theme(legend.title=element_text(size=6.5),
          legend.text=element_text(size=6),
          legend.key.size=unit(0.3, 'cm'))

  mu <- plyr::ddply(df, "method",
                    plyr::summarise,
                    grp.mean=mean(CV, na.rm=T))
  mu <- mu %>% mutate(Label=prettyNum(round(grp.mean, 4)))

  # peptide-centric histogram
  p2 <- ggplot(data=df, aes(x=CV)) +
    geom_histogram(color="black", fill="white") +
    facet_grid(method~.) +
    geom_vline(data=mu,
               aes(xintercept=grp.mean,
                   color="red"),
               linetype="dashed") +
    geom_text(
      data=mu,
      aes(x=grp.mean,
          y=nrow(df)/length(unique(df$method))/5,
          label=paste0("\u03bc=", Label)
          ),
      size=3,
      hjust=-.1) +
    theme_bw() +
    ylab("Frequency") +
    theme(legend.position="none")

  g1 <- gridExtra::grid.arrange(p2,p1, ncol=2)

  # protein group plots begin
  df <- rbind(cv5,cv6)
  df$method <- factor(df$method, levels=c("Protein Group", "Batch Adjusted"))

  df[,1] <- as.numeric(df[,1])
  df[,2] <- as.numeric(df[,2])

  df <- df[-which(is.na(df[,1])),]
  df <- df[-which(is.na(df[,2])),]

  # protein group density plot
  p3 <- ggplot(data=df,aes(`Log2 Median`,CV)) +
    geom_point(shape=16, size=0.1, color="grey66", show.legend=FALSE) +
    stat_density2d(aes(fill=..level..,alpha=..level..),
                   geom='polygon', colour='black') +
    scale_fill_continuous(low="green", high="red")+
    geom_smooth(method="gam", formula=y~s(x, k=2),
                linetype=2, colour="red",se=F) +
    guides(alpha="none") +
    facet_grid(method~., scales="free_x") +
    commonTheme +
    theme(legend.title=element_text(size=6.5),
          legend.text=element_text(size=6),
          legend.key.size=unit(0.3, 'cm'))

  mu <- plyr::ddply(df, "method", plyr::summarise,
                    grp.mean=mean(CV, na.rm=T))
  mu <- mu %>% mutate(Label=prettyNum(round(grp.mean, 4)))

  # protein group histogram
  p4 <- ggplot(data=df, aes(x=CV)) +
    geom_histogram(color="black", fill="white") +
    facet_grid(method ~ .) + geom_vline(data=mu,
                                        aes(xintercept=grp.mean, color="red"),
                                        linetype="dashed") +
    geom_text(
      data=mu,
      aes(
        x=grp.mean,
        y=nrow(df)/length(unique(df$method))/5,
        label=paste0("\u03bc=", Label)
      ),
      size=3,
      hjust=-.1
    ) +
    theme_bw() + ylab("Frequency") +
    theme(legend.position="none")

  g2 <- gridExtra::grid.arrange(p4,p3, ncol=2)

  plot_grid(g1,g2, ncol=1, labels=c('A' ,'B'), scale=0.975)

  ggsave("fig/figure3.pdf", dpi=600, dev='pdf', height=10, width=15, units="in")

  dev.off()

  # warnings back on
  options(warn=w)

}



