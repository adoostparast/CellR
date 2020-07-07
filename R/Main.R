#' CellR
#' @import dplyr Seurat  edgeR  glmnet  Matrix MASS Metrics gtools text2vec
#' @param Bulk
#' @param SingleCell
#' @param Bulk1
#' @param MinCell
#' @param Mingene
#' @param Dimension
#' @param Alpha
#' @param Alamat
#' @param maxi
#' @param mini
#' @param Data
#' @param Proportion
#' @param Cell_type

#library(devtools)
#devtools::install_deps(dependencies = TRUE)

#' @export
Deconvolution<-function(Bulk, SingleCell, Bulk1, Alamat,MinCell,Mingene,Dimension,Alpha, maxi, mini){


  library(dplyr)
  library(Seurat)
  library(edgeR)
  library(glmnet)
  library(Matrix)


  pbmc.data<-SingleCell
  #pbmc.data<-read.table("FCortex.txt",header=TRUE)
  #Bulk<-read.table("Bulk.txt",header=TRUE)
  #Bulk1<-read.table("GTExExpressionSymbol.txt",header=TRUE)
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.genes = Mingene, project = "CellR")
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  pbmc <- AddMetaData(object = pbmc, metadata = Alamat, col.name=Cell)
  par(mfrow = c(1, 2))
  #GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
  #GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
  pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(mini, -Inf), high.thresholds = c(maxi, 0.05))
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)
  #pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
  #Image<-JackStrawPlot(object = pbmc, PCs = 1:16)


  #  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:Dimension, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  #  pbmc <- RunTSNE(object = pbmc, dims.use = 1:Dimension, do.fast = TRUE)

  pbmc <- SetAllIdent(object = pbmc, id = "Cell")
  Clusters<-pbmc@ident
  Levels<-levels(Clusters)
  Levels<-as.data.frame(Levels)
  Clusters<-as.data.frame(Clusters)
  Mahsool.Clusters<-Clusters

  cluster.markers<-matrix(0,1,6)
  cluster.markers<-as.data.frame(cluster.markers)
  colnames(cluster.markers)<-c("p_val", "avg_logFC","pct.1", "pct.2", "p_val_adj", "cluster" )
  cluster.markers<-cluster.markers[-c(1),]


  for (i in 1:nrow(Levels)){
    kk<-Levels[i,1]
    kk<-as.data.frame(kk)
    kkk<-kk[1,1]
    kkk<-as.character(kkk)
    print("Finding markers of the cluster",quote = FALSE)
    print(kkk,quote = FALSE)
    cluster1.markers <- FindMarkers(object = pbmc, ident.1 = Levels[i,1], min.pct = 0.25)
    a<-matrix(Levels[i,1],nrow(cluster1.markers),1)
    a<-as.data.frame(a)
    cluster1.markers$cluster<-a[,1]
    cluster.markers<-rbind(cluster.markers,cluster1.markers)
  }

  sss<-which(cluster.markers$avg_logFC>0)
  sss<-as.data.frame(sss)
  FinalMarkers<-cluster.markers[sss[,1],]
  wqwq<-which(rownames(FinalMarkers) %in% rownames(Bulk))
  wqwq<-as.data.frame(wqwq)
  FinalMarkers<-FinalMarkers[wqwq[,1],]


  #--------------------------------------scMM----------------------------------------------------

  cc<-pbmc@scale.data
  cc<-as.data.frame(cc)
  ssMM<-cc[rownames(FinalMarkers),]
  ssMM<-as.data.frame(ssMM)

  rownames(ssMM)<-rownames(FinalMarkers)
  colnames(ssMM)<-colnames(cc)

  #----------------------------------Bulk RNA-Seq ----------------------------------------------
  bulkTamiz<-cpm(Bulk, normalized.lib.sizes=TRUE, log=FALSE)
  bulkTamiz<-as.data.frame(bulkTamiz)
  bulkTamizFinal<-bulkTamiz[rownames(FinalMarkers),]
  bulkTamizFinal<-as.data.frame(bulkTamizFinal)
  rownames(bulkTamizFinal)<-rownames(FinalMarkers)

  GTExNormalized<-Bulk1
  #GTExNormalized<-as.data.frame(GTExNormalized)
  Mean<-rowMeans(Bulk1)
  Mean<-as.data.frame(Mean)
  Temp<-rep(Mean,nrow(GTExNormalized),ncol(GTExNormalized))
  Temp<-as.data.frame(Temp)

  Mian<-(GTExNormalized-Temp)
  Mian<-as.data.frame(Mian)
  rownames(Mian)<-rownames(GTExNormalized)
  colnames(Mian)<-colnames(GTExNormalized)

  STD<-matrix(0,nrow(GTExNormalized),1)
  STD<-as.data.frame(STD)

  for (i in 1:nrow(GTExNormalized)){
    STD[i,1]<-sd(Mian[i,], na.rm = TRUE)}

  Zero<-matrix(0.00001,nrow(STD),1)
  Zero<-as.data.frame(Zero)
  STD<-Zero+STD

  Jam<-0

  for (k in 1:nrow(STD)){
    Jam<-(Jam+(1/STD[k,1]))}

  Final<-matrix(0,nrow(GTExNormalized),1)
  Final<-as.data.frame(Final)

  for (h in 1:nrow(STD)){
    Final[h,1]<-1+(1/(STD[h,1]*Jam))}

  rownames(Final)<-rownames(GTExNormalized)


  #------------------------------Optimization problem----------------------------------------
  ind<-matrix(0,nrow(ssMM),1)

  for (i in 1:nrow(ssMM)){
    if (is.na(ssMM[i,1])==TRUE){
      ind[i,1]<-1
    }
  }
  ind<-as.data.frame(ind)
  www<-which(ind[,c(1)]==0)
  www<-as.data.frame(www)
  FinalssMM<-ssMM[www[,1],]
  Finalbulk<-bulkTamizFinal[rownames(FinalssMM),]
  Finalbulk<-as.data.frame(Finalbulk)
  rownames(Finalbulk)<-rownames(FinalssMM)
  size<-ncol(Finalbulk)
  Finalbulki<-matrix(0,nrow(Finalbulk),1)
  Finalbulki<-as.data.frame(Finalbulki)
  rownames(Finalbulki)<-rownames(Finalbulk)

  GtexTemp<-Final[rownames(FinalssMM),]
  GtexTemp<-rep(GtexTemp,1,ncol(FinalssMM))
  #FinalssMM<-FinalssMM*GtexTemp
  #for (jiji in 1:ncol(FinalssMM)) {
  #  FinalssMM[,jiji]<-FinalssMM[,jiji]%*%GtexTemp}


  FinalssMM<-as.matrix(FinalssMM)
  #Finalbulk<-as.matrix(Finalbulk)
  UniqueClusters<-unique(Clusters)
  row.names(UniqueClusters)<-UniqueClusters[,1]
  OutputFinal<-matrix(0,nrow(UniqueClusters),ncol(Finalbulk))
  OutputFinal<-as.data.frame(OutputFinal)
  rownames(OutputFinal)<-row.names(UniqueClusters)

  OutputFreq<-OutputFinal   #Jadid

  for (Xha in 1:ncol(Finalbulk)){

    Finalbulki<-as.data.frame(Finalbulki)
    Finalbulki<-Finalbulk[,Xha]
    #  colnames(Finalbulki)<-colnames(Finalbulk[,Xha])
    Finalbulki<-as.matrix(Finalbulki)
    #-----------------------------------------------------------------------------------


    fit = glmnet(FinalssMM, Finalbulki, family = "gaussian",alpha = Alpha, nlambda = 100)
    CEO<-coef(fit,s=1)

    ttttt<-colnames(Finalbulki)
    ttttt<-as.data.frame(ttttt)

    Covariates<-matrix(0,nrow(Clusters), ncol(Finalbulki))
    Covariates<-as.data.frame(Covariates)
    rownames(Covariates)<-rownames(Clusters)
    colnames(Covariates)<-colnames(Finalbulki)


    tmp<-as.vector(CEO)
    tmp <- tmp[2:length(tmp)]
    Covariates[,1] <- tmp


    UniqueClusters<-unique(Clusters)
    row.names(UniqueClusters)<-UniqueClusters[,1]
    Output<-matrix(0,nrow(UniqueClusters),ncol(Covariates))
    Output<-as.data.frame(Output)
    row.names(Output)<-rownames(UniqueClusters)
    colnames(Output)<-colnames(Covariates)

    Frequency_Matrix<-Output  # Jadid

    for (i in 1:ncol(Covariates)){
      Tempo<-Covariates[,i]
      Tempo<-as.data.frame(Tempo)
      rownames(Tempo)<-rownames(Covariates)
      Tempo<-cbind(Tempo,rownames(Covariates))
      TemInd<-which(Tempo[,1]>0)
      Tempo<-Tempo[c(TemInd),]
      Tempo[,2]<-Clusters[rownames(Tempo),1]
      oo<-as.data.frame(table(Tempo[,2]))
      rownames(oo)<-oo$Var1
      Percent<-oo$Freq/sum(oo$Freq)*100
      oo$Var1<-Percent
      Output[,i]<-oo[rownames(Output),1]

      Frequency_Matrix[,i]<-oo[rownames(Output),c("Freq")] # Jadid


    }

    OutputFinal[,Xha]<-Output
    OutputFreq[,Xha]<-Frequency_Matrix # Jadid

    }

  colnames(OutputFinal)<-colnames(Bulk)
  colnames(OutputFreq)<-colnames(Bulk)  # Jadid

  FinalOutput<-list("Proportion"=OutputFinal,"Clusters"=Clusters, "Frequency"=OutputFreq)    # Taghir
  return(list("Proportion"=OutputFinal, "Markers"= FinalMarkers,"Clusters"=Clusters, "Frequency"=OutputFreq))   # Taghir
  #  return(list("Proportion"=OutputFinal))
}




Expression_estimate<-function(Data, Proportion, Frequency, Reference, Cell_type){  # Jadid

  library(MASS)
  library(Metrics)
  library(gtools)
  library(text2vec)
  library(edgeR)
  #Lambda=.5

  #Data<-cpm(Bulk, normalized.lib.sizes=TRUE, log=FALSE)
  Frequency<-t(Frequency)
  aa<-do.call("rbind", replicate(nrow(Data), colSums(Frequency), simplify = FALSE))  # Jadid
  Data<-Data/aa
  Data<-Data[,rownames(Proportion)]
  Reference<-Reference[rownames(Data),]

  T_max=300
  T_Min=1
  Temp=T_max
  Alpha=0.05
  Intermediate<-c()
  Counter=1

  Exp<-matrix(0,nrow(Proportion),nrow(Data))
  Exp<-as.data.frame(Exp)
  rownames(Exp)<-rownames(Proportion)  # Jadid
  colnames(Exp)<-rownames(Data) # Jadid


  for (i in 1:nrow(Data)){
    cat("\nGene =", i)

    #X<-cbind(matrix(1,nrow(Proportion),1),Proportion)
    #Y<-t(Data[c(i),])
    #XX<-as.matrix(X)

    Mean=Reference[c(i),]     # Mean_Ast >-Mean

    theta = runif(1, 0, 1)                  # theta_Ast >-theta


    Matris=matrix(0,ncol(Proportion),2)
    Matris<-as.data.frame(Matris)
    Matris[,c(1)]<-t(Mean)
    Matris[,c(2)]<-theta
    Matris_new<-Matris

    Ast<-as.data.frame(matrix(0,nrow(Proportion),ncol(Proportion)))
    Ast_New<-Ast

    for (w in 1:nrow(Matris)){
      Ast[,c(w)]<-rnegbin(ncol(Data), mu = Matris[c(w),c(1)], theta = Matris[c(w),c(2)])

    }

    # Ast<-as.data.frame(rnegbin(ncol(Data), mu = Mean, theta = theta))

    Intermediate<-Ast
    # colnames(Intermediate)<-Cell_type
    Jam<-Proportion*Intermediate
    ee<-as.data.frame(rowSums(Jam))
    Error_old<-rmse(as.numeric(as.matrix(ee)),as.numeric(as.matrix(Data[c(i),]))) #+Lambda*sum(Matris[,c(1)]^2)



    while (Temp>T_Min){



      for (k in 1:200){

        Prob=sample(1:2, 1)

        if (Prob==1){


          Matris_new[,c(1)]<-Matris[,c(1)]+runif(nrow(Matris), -.5, .5)*Matris[,c(1)]
          Matris_new[,c(2)]<-Matris[,c(2)]+runif(nrow(Matris), 0, 1)*Matris[,c(2)]

          # Mean_New=Mean+runif(1, -.5, .5)*Mean
          # theta_New = theta+runif(1, -.3, .3)*theta
        }

        if (Prob==2) {

          Matris_new[,c(1)]<-as.data.frame(runif(nrow(Matris), mean(as.numeric(Data[c(i),]))-sd(as.numeric(Data[c(i),])), mean(as.numeric(Data[c(i),]))+sd(as.numeric(Data[c(i),]))))
          Matris_new[,c(2)]<-as.data.frame(runif(nrow(Matris),0,1))
          #Matris_new[,c(2)]<-as.data.frame(sample(1:5, nrow(Matris),replace=TRUE))

          # Mean_New=runif(1, mean(as.numeric(Data[c(i),]))-sd(as.numeric(Data[c(i),])), mean(as.numeric(Data[c(i),]))-sd(as.numeric(Data[c(i),])))
          # theta_New = sample(1:5, 1)
        }

        #        if (Prob==3) {

        #        rows <- sample(nrow(Matris))
        #         Matris_new[,c(1)] <- Matris[rows, c(1)]
        #         Matris_new[,c(2)] <- Matris[rows, c(2)]


        #      }

        for (w in 1:nrow(Matris)){
          Ast_New[,c(w)]<-rnegbin(ncol(Data), mu = Matris_new[c(w),c(1)], theta = Matris_new[c(w),c(2)])

        }

        # Ast_New<-as.data.frame(rnegbin(ncol(Data), mu = Mean_New, theta = theta_New))

        Intermediate1<-Ast_New
        #colnames(Intermediate1)<-Cell_type
        Jam1<-Proportion*Intermediate1
        ee1<-as.data.frame(rowSums(Jam1))
        Error_New<-rmse(as.numeric(as.matrix(ee1)),as.numeric(as.matrix(Data[c(i),]))) # +Lambda*sum(Matris_new[,c(1)]^2)

        if (Error_New<Error_old){


          Matris<-Matris_new
          #Mean<-Mean_New
          #theta<-theta_New

          #Matris[c(1),c(1)]<-Mean
          #Matris[c(1),c(2)]<-theta


          cat("\nTemp =", Temp)
          cat("\nK =", k)
          cat("\nError =", Error_New, "\n")

          Error_old<-Error_New}

        else {

          P=1-exp(Temp*(Error_old-Error_New))
          if(P<0.01){

            Matris<-Matris_new
            # Mean<-Mean_New
            #theta<-theta_New

            # Matris[c(1),c(1)]<-Mean
            # Matris[c(1),c(2)]<-theta

            Error_old<-Error_New

          }
        }


      }

      Counter=Counter+1
      Temp=T_max*(1-Alpha)^Counter

    }
    Counter=1
    Temp=T_max



    Re<-Matris[c(Cell_type),1]
    Varoon<-rep(Re, each=nrow(Proportion))
    Out<-Frequency[,c(Cell_type)]*Varoon

    # Varoon<-t(as.data.frame(t(Matris)[,c(Cell_type)]))
    #Varoon<-Varoon[rep(seq_len(nrow(Varoon)), each = nrow(Proportion)), ]
    #Out<-Proportion[,c(Cell_type)]*Varoon
    #Out<-Varoon # Jadid

    Exp[,c(i)]<-Out
  }
  #FinalExp<-Exp[c(1),]
  #colnames(FinalExp)<-rownames(Data)
  return(Exp)
}

