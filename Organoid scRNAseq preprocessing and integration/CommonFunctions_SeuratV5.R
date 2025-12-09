#This is a collection of functions useful for Seurat objects.
#To load it run:
#source("commonFunctions.R")

#This is part of VlnPlot.median function
median.stat <- function(x){
	out <- quantile(x, probs = c(0.5))
	names(out) <- c("ymed")
	return(out) 
	}

VlnPlot.median <- function(object, features, colour = "black", pt.size = 0, ncol = NULL, legend = FALSE, vln.colors=NULL) {
	#Tris function creates VlnPlot plus a median 
	#VlnPlot.median(object, c("IL1B","CLEC10A"))
	
	myplots <- vector("list")
	
	#Create a plot for each gene and add median
	for (gene in features) {
		#print(gene)
		if (legend == TRUE) {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size, cols = vln.colors) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour)
		} else {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size, cols = vln.colors) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour) + NoLegend()
		}
	}
	#patchwork function to combine multiple plots
	patchwork::wrap_plots(myplots, ncol=ncol)
	}


SoupX.clean.from.CellRanger <- function(cellranger.folder) {
  #Clean your scRNAseq data from ambient contamination
  #Use:
  #clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")
  
  # if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
  #   stop("filtered_feature_bc_matrix folder didn't find")
  # }
  # 
  # if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
  #   stop("raw_feature_bc_matrix folder didn't find")
  # }
  # 
  # if (dir.exists(paste(cellranger.folder,"analysis",sep = "/"))==FALSE){
  #   stop("analysis folder didn't find")
  # }
  
  #Load data
  sc = SoupX::load10X(cellranger.folder)
  #Estimante the contamination
  sc = SoupX::autoEstCont(sc)
  
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  
  #Filterout the contamination
  out = SoupX::adjustCounts(sc)
  return(out)
}

SoupX.on.Seurat <- function(Seurat.object, cellranger.folder, min.cells = 5, min.features = 50) {
  #Clean your scRNAseq data from ambient contamination
  #Seurat object should be UNFILTERED (MT genes and doublets), and seurat_clusters should be present
  #Use:
  #clean.data = SoupX.clean(Seurat.object, cellranger.folder="/my/cellranger/folder/")
  
  if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
    stop("filtered_feature_bc_matrix folder didn't find")
  }
  
  if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
    stop("raw_feature_bc_matrix folder didn't find")
  }
  
  if (any(colnames(Seurat.object@meta.data) == "seurat_clusters")) {
    print("seurat_clusters present...")
  } else {
    stop("seurat_clusters not present")
  }
  
  raw.matrix = Seurat::Read10X(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))
  filt.matrix = Seurat::Read10X(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))
  
  sc  <- SoupX::SoupChannel(raw.matrix, filt.matrix)
  
  meta    <- Seurat.object@meta.data
  umap    <- Seurat.object@reductions$umap@cell.embeddings
  sc  <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- SoupX::setDR(sc, umap)
  head(meta)
  
  #With defined clusters, run the main SoupX function, calculating ambient RNA profile.
  sc  <- SoupX::autoEstCont(sc)
  
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  
  adj.matrix = SoupX::adjustCounts(sc, roundToInt = T)
  
  New.Seurat.object = Seurat::CreateSeuratObject(counts = adj.matrix, min.cells = 5, min.features = 50, meta.data = meta)
  return(New.Seurat.object)
}

make.add.meta <- function(Seurat.Object, metadata, return.only.table=FALSE, verbose=FALSE) {
  #From a metadata table and a Seurat.Object, 
  #it creates a proper metadata table with cell barcodeID and add it to Seurat.object 
  # Use:
  # Seurat.object = make.add.meta(Seurat.Object, metadata)
  
  if (verbose) {
    print(head(metadata))
  }
  
  #if there is only 1 row, add it to the whole Seurat.Object
  if (nrow(metadata)==1) {
    if (verbose) {
      print("nrow(metadata)==1")
    }
    df.cells <- data.frame(row.names = colnames(Seurat.Object))
    for (name in colnames(metadata)) {
      #print(name)
      df.cells[name]=metadata[[name]]
    }
  #if there is only 1 column it's a list  
  } else if (ncol(metadata)==1) {
    if (verbose){
      print(colnames(metadata))
      print(names(metadata))
    }
    
    #check if length of idents is different of lenght metadata
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata)))!=0 || length(setdiff(rownames(metadata),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata rows are not matching.")
    }
    
    df.cells <- data.frame()
    
    #Select the Idents matching “Name” in metadata
    #Idents(Seurat.Object) <- "Condition_name"
    
    #For each idents within the seurat object (cluster or sample you want)
    for (name in unique(Idents(Seurat.Object)))
    {
      print(name)
      
      #Select the cells and put them in another df
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))
      
      #Select the row of interest from metadata, corresponding to the metadata to add to those cells
      meta_row=metadata[rownames(metadata) == name,]
      new_df[colnames(metadata)]=meta_row
      df.cells=dplyr::bind_rows(df.cells, new_df)
      }
      #merge in a big df
  #if else is a dataframe
  } else {
    if (verbose) {
      print("nrow(metadata)!=1")
    }
    #check if length of idents is different of lenght metadata
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata)))!=0 || length(setdiff(rownames(metadata),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata rows are not matching.")
    }
    
    df.cells <- data.frame()
    
    #Select the Idents matching “Name” in metadata
    #Idents(Seurat.Object) <- "Condition_name"
    
    #For each idents within the seurat object (cluster or sample you want)
    for (name in unique(Idents(Seurat.Object)))
    {
      print(name)
      
      #Select the cells and put them in another df
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))
      
      #Select the row of interest from metadata, corresponding to the metadata to add to those cells
      meta_row=metadata[rownames(metadata) == name,]
      if (verbose) {
        #print(metadata)
        #print(rownames(metadata))
        print("meta_row:")
        print(meta_row)
        #print(nrow(meta_row))
        print(colnames(meta_row))
        print(ncol(meta_row))
      }
      if (is.null(colnames(meta_row))) {
        stop("Contact Dom")
      }
      if (ncol(meta_row)==1) {
        stop("There is only 1 column...contact dom, need to be fixed")
      }
      for (col_name in colnames(meta_row)){
        #add the specific metadata you need
        if (verbose) {
          print("col_name is:")
          print(col_name)
        }
        new_df[col_name]=meta_row[[col_name]]
      }
      #merge in a big df
      df.cells=dplyr::bind_rows(df.cells, new_df)
    }
  }
  
  #Finally add to the object
  
  if (return.only.table==TRUE) {
    return(df.cells)
  } else {
    if (verbose) {
      print("Adding to Seurat...")
    }
    Seurat.Object <- AddMetaData(Seurat.Object, metadata = df.cells)
    return(Seurat.Object)
  }
}

Mark.cells <- function(Seurat.object, nFeature.low = 250, nFeature.high = 5000, mt.high = 25, nCount.high = 18000) {
  #Mark cells with low/high nFeature, nCount.high, mt genes,
  #Example use: Seurat.Object <- Mark.cells(Seurat.Object)
  
  # Known issues, to improve:
  # this function works on "percent.mt", "nFeature_RNA" and "nCount_RNA" metadata columns...need to be more flexible
  
  poscells.high.mt <- WhichCells(Seurat.object, expression = percent.mt > mt.high)
  Seurat.object$high.mt<- ifelse(colnames(Seurat.object) %in% poscells.high.mt, "Pos", "Neg")
  print("high mt cells:")
  print(table(Seurat.object$high.mt))
  
  poscells.high.ft <- WhichCells(Seurat.object, expression = nFeature_RNA > nFeature.high)
  Seurat.object$high.nFeature <- ifelse(colnames(Seurat.object) %in% poscells.high.ft, "Pos", "Neg")
  print("high nFeature cells:")
  print(table(Seurat.object$high.nFeature))
  
  poscells.low.ft <- WhichCells(Seurat.object, expression = nFeature_RNA < nFeature.low)
  Seurat.object$low.nFeature <- ifelse(colnames(Seurat.object) %in% poscells.low.ft, "Pos", "Neg")
  print("low nFeature cells:")
  print(table(Seurat.object$low.nFeature))
  
  poscells.high.nc <- WhichCells(Seurat.object, expression = nCount_RNA > nCount.high)
  Seurat.object$high.nCount <- ifelse(colnames(Seurat.object) %in% poscells.high.nc, "Pos", "Neg")
  print("high UMIs:")
  print(table(Seurat.object$high.nCount))
  
  return(Seurat.object)
}

find.significant.PCs <- function(Seurat.Object, variance=0.95, st.dev=0.05, reduction="pca") {
  ##Find a min PCA significant (90% variance & st.dev<5%) to use
  ##The Seurat.Object need "pca" slot filled
  ##Warning: it calculates 90% on n.pca calculated!! It will gives different results if you calculate 30 or 50 pca dims
  ##Example use: 
  ## min.pca = find.significant.PCs(Seurat.Object)
  
  ##Describe the code
  
  if (variance>=1) {
    simpleError("variance over 100%")
  }
  variance=variance*100
  st.dev=st.dev*100
  
  stdv <- Seurat.Object[[reduction]]@stdev
  sum.stdv <- sum(Seurat.Object[[reduction]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > variance & percent.stdv < st.dev)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 1-(variance/100)), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  print(min.pc)
  return(min.pc)
}

Calc.Perc.Features <- function(Seurat.object, mt.pattern = "^MT-", hb.pattern = "^HB[^(P)]", ribo.pattern = ("RPS|RPL"), MALAT1.name="MALAT1", plot.name="") {
  #v0.2
  Seurat.object@misc$cell.recovered = ncol(Seurat.object)
  Seurat.object[["percent.mt"]] <- PercentageFeatureSet(Seurat.object, pattern = mt.pattern)
  Seurat.object[["percent.hb"]] <- PercentageFeatureSet(Seurat.object, pattern = hb.pattern)
  Seurat.object[["percent.ribo"]] <- PercentageFeatureSet(Seurat.object, pattern = ribo.pattern)
  Seurat.object[["percent.MALAT1"]] <- PercentageFeatureSet(Seurat.object, pattern = MALAT1.name)
  plot1 <- FeatureScatter(Seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(Seurat.object, feature1 = "percent.MALAT1", feature2 = "percent.mt")
  plot4 <- VlnPlot(Seurat.object, features = c("percent.ribo", "percent.hb"), ncol = 2)
  plot <- ((plot1 + plot2) / (plot3 + plot4))
  plot = plot + plot_annotation(title = plot.name, theme = theme(plot.title = element_text(hjust = 0.5)))
  print(plot)
  return(Seurat.object)
}

QC.n.mad <- function(Seurat.object, n.mad=4) {
  #Based on https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html code
  Cell.QC.Stat <- Seurat.object@meta.data
  print(nrow(Cell.QC.Stat))
  # high and low median absolute deviation (mad) thresholds to exclude outlier cells
  max.mito.thr <- median(Cell.QC.Stat$percent.mt) + n.mad*mad(Cell.QC.Stat$percent.mt)
  min.mito.thr <- median(Cell.QC.Stat$percent.mt) - n.mad*mad(Cell.QC.Stat$percent.mt)
  
  #Plot
  p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[2])," cells removed\n", as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[1])," cells remain"), x = 6000, y = 0.1)
  
  Cell.QC.Stat <- Cell.QC.Stat %>% dplyr::filter(percent.mt < max.mito.thr) %>% dplyr::filter(percent.mt > min.mito.thr)
  
  # Set low and hight thresholds on the number of detected genes
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  # Set hight threshold on the number of transcripts
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nCount_RNA))
  
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)
  
  p3 <- p1 / p2
  #p1=ggExtra::ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)
  #p2=ggExtra::ggMarginal(p2, type = "histogram", fill="lightgrey")
  #print(typeof(p1))
  #p3 <- p1 / p2
  p3 <- p3 + plot_annotation(title = names(Seurat.object), theme = theme(plot.title = element_text(hjust = 0.5)))
  print(p3)
  
  # Filter cells based on these thresholds
  
  Cell.QC.Stat <- Cell.QC.Stat %>% dplyr::filter(log10(nFeature_RNA) > min.Genes.thr) %>% dplyr::filter(log10(nFeature_RNA) < max.Genes.thr) %>% dplyr::filter(log10(nCount_RNA) < max.nUMI.thr)
  print(nrow(Cell.QC.Stat))
  print("########")
  Seurat.object <- subset(Seurat.object, cells = rownames(Cell.QC.Stat))
  return(Seurat.object)
}