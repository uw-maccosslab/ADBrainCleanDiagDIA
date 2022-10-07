#' Sci Data Paper Supplementary Figure 5
#' 
#' @import lme4
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics axis barplot par text
#' @return supplfig5.pdf
#' @export

get_supfig5 <- function() {
  
  printPVCA = function(theDataMatrix, metaDT, main){
    
    pct_threshold = .4 # Amount of variability desired to be explained by the principal components.
    
    theDataMatrix[is.na(theDataMatrix)] = 0
    dataRowN <- nrow(theDataMatrix)
    dataColN <- ncol(theDataMatrix)
    
    theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)
    exp_design <- metaDT
    expDesignRowN <- nrow(exp_design)
    expDesignColN <- ncol(exp_design)
    myColNames <- names(exp_design)
    
    theDataCor <- cor(theDataMatrixCentered)
    
    eigenData <- eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n <- length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues /eigenValuesSum 
    
    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
      my_sum_2  = my_sum_2 - percents_PCs[i]
      if ((my_sum_2) <= pct_threshold ){
        my_counter_2 = my_counter_2 + 1
      }
      
    }
    if (my_counter_2 < 3){
      pc_n  = 3
      
    }else {
      pc_n = my_counter_2 
    }
    
    pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
      for (j in 1:expDesignRowN){
        mycounter <- mycounter + 1
        pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        
      }
    }
    
    AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
    
    Data <- cbind(AAA,pc_data_matrix)
    
    Data$Batch <- as.factor(Data$Batch)
    Data$Condition <- as.factor(Data$Condition)
    Data$Age <- as.numeric(Data$Age)
    Data$Sex <- as.factor(Data$Sex)
    
    op <- options(warn = (-1)) 
    effects_n = (expDesignColN - 2) + ((expDesignColN - 2)*(((expDesignColN - 2)-1)))/2 + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    
    for (i in 1:pc_n){
      y = (((i-1)*expDesignRowN)+1)
      
      Rm1ML <- lmer(pc_data_matrix ~ (1|Batch) + (1|Condition) + (1|Age) + (1|Sex), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit)
      
      randomEffects = as.data.frame(VarCorr(Rm1ML))
      for (j in 1:effects_n){
        randomEffectsMatrix[i,j] = as.numeric(randomEffects[j,4])
      }
      
    }
    
    effectsNames <- randomEffects[,1]
    
    randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
      mySum = sum(randomEffectsMatrix[i,], na.rm = T)
      for (j in 1:effects_n){
        randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
      }
    }
    
    randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
      weight = eigenValues[i]/eigenValuesSum
      for (j in 1:effects_n){
        randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
      }
    }
    
    randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums, na.rm = T)
    randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
    
    for (j in 1:effects_n){
      randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
      
    }
    

    bp <- barplot(randomEffectsMatrixWtAveProp[!is.na(randomEffectsMatrixWtAveProp)],  
                  xlab = "Effects", ylab = "Weighted average proportion variance", 
                  main = main,
                  ylim= c(0,1.1),col = c("blue"), las=2)
    
    axis(1, at = bp[1:length(effectsNames)], labels = effectsNames, xlab = "Effects", cex.axis = 1, las=1)
    values = randomEffectsMatrixWtAveProp
    new_values = round(values , 3)
    text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) # place numbers on top of bars 

  }
  
  
  # warnings off
  w <- getOption("warn")
  options(warn=-1)
  
  if( !exists("fig"))
    dir.create("fig")
  
  # load data
  data(protGrp, envir=environment())
  data(noBatchProt, envir=environment())
  data(metaDT, envir=environment())
  
  protGrp <- as.matrix(protGrp)[,order(colnames(protGrp))]
  noBatchProt <- as.matrix(noBatchProt)[,order(colnames(noBatchProt))]
  metaDT = metaDT[order(metaDT$Replicate),]
  protGrp[protGrp <= 0] <- NA
  noBatchProt[noBatchProt <= 0] <- NA
  
  protGrp = protGrp[complete.cases(protGrp),]
  noBatchProt = noBatchProt[complete.cases(noBatchProt),]

  filename = "supplfig5"
  
  pdf(
    paste0("fig/",filename,".pdf"),
    height = 9,
    width = 15,
    compress = F
  )
  
  par(mfrow = c(1,2))
  
  ## Pre-Batch Adjust PVCA #####################################################
  printPVCA(protGrp, metaDT, "Pre-Adjustment")
  
    
  ## Post Batch Adjust PVCA #####################################################
  printPVCA(noBatchProt, metaDT, "Post-Adjustment")
  
  
  dev.off()
  
  # warnings back on
  options(warn = w)
}
