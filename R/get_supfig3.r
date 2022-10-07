#' Sci Data Paper Supplementary Figure 3
#'
#' @import mdatools
#'
#' @return supplfig3.pdf
#' @export
get_supfig3 = function() {

  # warnings off
  w <- getOption("warn")
  options(warn = -1)

  if( !exists("fig"))
    dir.create("fig")

  # load data
  data(noBatchPep, envir = environment())
  data(metaDT, envir = environment())
  
  # rename condition groups
  metaDT$Condition = gsub("Path", "ADNC", gsub("Control ", "HCF - ", 
                            gsub("AD", "ADD", metaDT$Condition)))
  dat = noBatchPep

  # wrangle
  dat[dat == 0] = NA
  dat[is.na(dat)] = 0
  dat = dat[complete.cases(dat),]
  dat = as.matrix(dat)[, match(metaDT$Replicate, colnames(dat))]
  metaDT$Age[is.na(metaDT$Age)] = 0
  metaDT$Condition = gsub(
    "SMTG Reference",
    "Interbatch QC",
    gsub("Brain Reference",
         "Interexperiment QC",
         metaDT$Condition)
  )

  m <- mdatools::pca(t(as.data.frame(dat)), center = T, scale = T)

  pdf(
    "fig/supplfig3.pdf",
    height = 9,
    width = 15,
    compress = F
  )

  par(mfrow = c(2, 2))
  plotScores(m,
             show.labels = TRUE,
             main = "Replicate",
             ylim = c(-200, 210))

  p = plotScores(
    m,
    c(1, 2),
    cgroup = metaDT$Age,
    main = "Age",
    ylim = c(-200, 210)
  )

  # scores plot colored by the factor created above and confidence ellipses
  p = plotScores(
    m,
    c(1, 2),
    cgroup = metaDT$Batch,
    main = "Batch",
    ylim = c(-200, 210)
  )
  plotConfidenceEllipse(p)

  new_colors = cbind(
    Condition = c(
      "AutoDom ADD",
      "Sporadic ADD",
      "HCF - High ADNC",
      "HCF - Low ADNC",
      "Interexperiment QC",
      "Interbatch QC"
    ),
    Colors = c(
      "#ff4500",
      "#ffa500",
      "#000080",
      "#1874cd",
      "#8c8c8c",
      "#cccccc"
    )
  )

  p = plotScores(
    m,
    c(1, 2),
    cgroup = factor(
      metaDT$Condition,
      levels = c(
        "AutoDom ADD",
        "Sporadic ADD",
        "HCF - High ADNC",
        "HCF - Low ADNC",
        "Interexperiment QC",
        "Interbatch QC"
      )
    ),
    main = "Condition",
    colmap = new_colors[, "Colors"],
    ylim = c(-200, 210)
  )
  plotConfidenceEllipse(p)

  dev.off()

  # warnings back on
  options(warn = w)
}


