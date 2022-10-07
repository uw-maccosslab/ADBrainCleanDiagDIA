#' Sci Data Paper Supplementary Figures 1 and 2
#'
#' @import readr
#' @import ggplot2
#' @import dplyr
#' @import ggridges
#'
#' @return supplfigure1.pdf;supplFigure2.pdf
#' @export

get_supfig1_2=function() {
  
  # warnings off
  w <- getOption("warn")
  options(warn=-1)
  
  if( !exists("fig"))
    dir.create("fig")

  files=list.files(path = "inst/extdata/",pattern = paste0("Pivot.csv","$") ,full.names = TRUE)
  raw.batch = lapply(files, function(x) read_csv(x))
  
  for (i in 1:length(raw.batch)){
    colnames(raw.batch[[i]]) = sapply(strsplit(colnames(raw.batch[[i]]), " "), "[", 1)
    colnames(raw.batch[[i]]) = sapply(strsplit(colnames(raw.batch[[i]]), "-"), "[", 1)
    if (length(grep("Peptide", colnames(raw.batch[[i]]))) > 1)
      raw.batch[[i]] = subset(raw.batch[[i]], select = -Peptide)
    if (length(grep("PeptideModifiedSequence", colnames(raw.batch[[i]]))) > 0)
      colnames(raw.batch[[i]])[which(colnames(raw.batch[[i]]) == "PeptideModifiedSequence")] = "Peptide"
    if (length(grep("ProteinName", colnames(raw.batch[[i]]))) > 0)
      colnames(raw.batch[[i]])[which(colnames(raw.batch[[i]]) == "ProteinName")] = "Protein"
    uniquePeptide = paste0(raw.batch[[i]]$Peptide, "@", raw.batch[[i]]$Protein)
    raw.batch[[i]] = raw.batch[[i]][,-c(1:2)]
    rownames(raw.batch[[i]]) = uniquePeptide
  }
  
  
  raw.batch = lapply(raw.batch, function(x) {x$uniquePeptide = rownames(x); return(x)})
  combo = full_join(raw.batch[[1]],raw.batch[[2]], by = "uniquePeptide")
  if (length(raw.batch) > 2)
    for (i in 3:length(raw.batch))
      combo = full_join(combo,raw.batch[[i]], by = "uniquePeptide")
  
  
  rnames = combo$uniquePeptide
  df = combo[, -which(colnames(combo) == "uniquePeptide")]
  df = df[,order(colnames(df))]
  
  df[df==0] = NA
  df[log2(df) < -20] = NA
  df = log2(df)
  
  rownames(df) = rnames
  df = as.matrix(df)
  
  
  metaFiles=list.files(path = "inst/extdata", pattern = "meta.csv$", full.names = TRUE)
  meta.dt = lapply(metaFiles, function(x) read_csv(x))
  
  for (i in 1:length(meta.dt)){
    meta.dt[[i]]$Region = sapply(strsplit(metaFiles, "-"), "[", 3)[i]
    meta.dt[[i]]$Batch  = paste0(meta.dt[[i]]$Region, meta.dt[[i]]$Batch)
  }
  
  temp = data.frame()
  for (i in 1:length(meta.dt))
    temp = rbind(temp, meta.dt[[i]])
  meta.dt = temp
  
  meta.dt$Replicate = meta.dt$`Sample Label`
  
  meta.dt = meta.dt[colnames(meta.dt) != "Sample Label"]
  
  # rename condition groups
  meta.dt$Condition = gsub("Path", "ADNC", gsub("Control ", "HCF - ", 
                            gsub("AD", "ADD", meta.dt$Condition)))
  
  df = df[,match(meta.dt$Replicate, colnames(df))]
  
  
  df.long = reshape2::melt(df, value.name = "Abundance")
  names(df.long) = c("Peptide", "Replicate", "Abundance")
  df.long = merge(df.long, meta.dt[,c("Replicate", "Condition", "Batch", "Region")], by="Replicate")
  
  df.long$Batch = gsub("Caudate|Hipp|IPL|SMTG","", df.long$Batch)
  
  df.long$Region = factor(df.long$Region, levels = c("SMTG","Hipp","IPL","Caudate"))
  
  p1 = ggplot(df.long, aes(x = Replicate, y = Abundance, fill = NULL)) +
    geom_boxplot() + theme(legend.position = "none") + facet_wrap(~Region, scales = "free_x") +
    labs(y="Log2 Abundance") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = rel(0.7)))
  
  p2 = ggplot(df.long, aes(x = Abundance, y = Batch, fill = NULL)) +
    geom_density_ridges(scale = 0.85,
                        jittered_points = F,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '|', point_size = 0.9, point_alpha = 0.1, alpha = 0.7,
    ) + facet_grid(~Region) + labs(x="Log2 Abundance") + theme(legend.position = "none")
  
  g1 <- gridExtra::grid.arrange(p1,p2, ncol=1)
  
  ggsave("fig/supplfig1.pdf",g1, dpi=600, dev='pdf',
         height=9, width=15, units="in")
  dev.off()
  
  
  df.norm = df
  for (i in 1:length(unique(meta.dt$Region))) {
    idx = which(match(meta.dt$Region, unique(meta.dt$Region)) == i)
    df.norm[,idx] = normalize_median(df.norm[,idx])
  }
  
  df.long = reshape2::melt(df.norm, value.name = "Abundance")
  names(df.long) = c("Peptide", "Replicate", "Abundance")
  df.long = merge(df.long, meta.dt[,c("Replicate", "Condition", "Batch", "Region")], by="Replicate")
  
  df.long$Batch = gsub("Caudate|Hipp|IPL|SMTG","", df.long$Batch)
  
  df.long$Region = factor(df.long$Region, levels = c("SMTG","Hipp","IPL","Caudate"))
  
  p3 = ggplot(df.long, aes(x = Replicate, y = Abundance, fill = NULL)) +
    geom_boxplot() + theme(legend.position = "none") + facet_wrap(~Region, scales = "free_x") +
    labs(y="Log2 Abundance") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = rel(0.7)))
  
  p4 = ggplot(df.long, aes(x = Abundance, y = Batch, fill = NULL)) +
    geom_density_ridges(scale = 0.85,
                        jittered_points = F,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '|', point_size = 0.9, point_alpha = 0.1, alpha = 0.7,
    ) + facet_grid(~Region) + labs(x="Log2 Abundance") + theme(legend.position = "none")
  
  g2 = gridExtra::grid.arrange(p3,p4, ncol=1)
  
  ggsave("fig/supplfig2.pdf",g2, dpi=600, dev='pdf',
         height=9, width=15, units="in")
  dev.off()
  
  # warnings back on
  options(warn=w)
}
