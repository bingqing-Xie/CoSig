#' Setting multi-threading
#' @param cores Number of cores to use, not to overload your computer
initializeMultiCores = function(cores=10){
	cl <- makeCluster(cores) 
	registerDoParallel(cl)
}

#' Compute co-expression networks given seurat object for each sample and
#' Perform unsupervised HOSVD on the combined 3D tensor co-expression network
#' @param prefix folder for raw output (automatic create raw folder)
#' @param sampleid column name in seurat containing sample info
#' @param recompute_glasso recompute the graphical lasso co-expression network
#' @param recompute_HOSVD recompute HOSVD for the co-expression network
#' @param ReComputeVarGene recompute the variable genes
#' @param mean_thres minimal threshold for mean expression of variable genes
#' @param return Decomposed tensor matrices , list of sample names , and list of gene names
#' @export
computeCoExp  =  function(obj, prefix = ".", sampleid = "Sample_ID",recompute_glasso = TRUE, recompute_HOSVD = TRUE, ReComputeVarGene = FALSE, mean_thres = 0.1){
	if(!dir.exists(file.path(prefix))){
		prefix = "."
	}
	dir.create(file.path(prefix, "raw"), showWarnings  =  FALSE)
	prefix = paste0(prefix,"/raw/")
	if(recompute_glasso){
		if(ReComputeVarGene){
			obj = FindVariableFeatures(obj,nfeatures  =  2000)
			hvf = HVFInfo(obj)[VariableFeatures(obj),]
			hvf = hvf[order(hvf$mean,decreasing = T),]
			varFeatures = intersect(rownames(subset(HVFInfo(obj),mean>mean_thres)),VariableFeatures(obj))
			write.csv(varFeatures,file = paste0(prefix,"varFeatures.csv"))
		}else{
			varFeatures = read.csv(file = paste0(prefix,"varFeatures.csv"),row.names = 1)$x
		}	
	
		objlist = SplitObject(obj,sampleid)
		sampleNames = names(objlist)
		write.csv(sampleNames,file = paste0(prefix,"sampleNames.csv"))
		nFeatures = length(varFeatures)
		nsamples = length(sampleNames)
		t_coexp = rand_tensor(modes = c(nFeatures,nFeatures,nsamples))
		coexp_mats = foreach(i = 1:length(objlist)) %dopar% {
			print(i)
			s = var(t(as.matrix(GetAssayData(objlist[[i]],slot = "counts")[varFeatures,] )))
			#coexp = glassoFast(s, rho = 0.1)	
			coexp = glasso(s, rho = 0.1)	
			wi  =  coexp$wi
			for(j in 1:nFeatures){
				wi[j,j:nFeatures] = wi[j:nFeatures,j]
			}
			return(wi)
		}

		for(i in 1:length(objlist)){
			t_coexp[,,i] = coexp_mats[[i]]
		}
		saveRDS(t_coexp,file = paste0(prefix,"coexp_networks_tensor.rds"))
	}else{
		varFeatures = read.csv(file = paste0(prefix,"_varFeatures.csv"),row.names = 1)[,1]
		sampleNames = read.csv(file = paste0(prefix,"_sampleNames.csv"),row.names = 1)[,1]		
		nFeatures = length(varFeatures)
		nsamples = length(sampleNames)
		t_coexp = readRDS(file = paste0(prefix,"coexp_networks_tensor.rds"))
	}
	
	if(recompute_HOSVD){
		cp_decomp <- cp(t_coexp, num_components  =  nsamples, max_iter  =  100)
		str(cp_decomp$U)
		saveRDS(cp_decomp,file = paste0(prefix,"HOSVD_cp_decomp.rds"))
	}else{
		cp_decomp = readRDS(file = paste0(prefix,"HOSVD_cp_decomp.rds"))
	}
	return(list(cp_decomp,t_coexp,sampleNames,varFeatures))
}

#' Extract the tensor component (TC) that correlated with the feature of interest,
#' and construct the network
#' @param prefix working directory prefix
#' @param cp_decomp decomposed matrices from HOSVD
#' @param t_coexp 3D tensor matrix for stacked co-expression networks for all samples
#' @param sampleNames names of the samples
#' @param varFeatures gene names
#' @param meta metadata for samples
#' @param featureOfInterest column name of the feature of interest in the metadata assuming the features are levels or numerical values
#' @param featureOfCategory column name of the grouping in the metadata: assuming the features are levels
#' @param cor_method correlation method
#' @param zthres threshold for z-score of loading values to identify significant driving genes per TC
#' @export
featureAnalysis =  function(prefix, cp_decomp,t_coexp,sampleNames,varFeatures, meta, featureOfInterest, featureOfCategory, cor_method  =  "spearman", zthres = 2){
	if(!dir.exists(file.path(prefix))){
		prefix = "."
	}
	dir.create(file.path(prefix, "result"), showWarnings  =  FALSE)
	prefix = paste0(prefix,"/result/")
	### Extract the sample x TC matrix and concatenate it with the metadata
	df2 = data.frame(cp_decomp$U[[3]])
	colnames(df2) = paste0("Comp",1:ncol(df2))
	df2$Sample_ID = sampleNames
	df2 = cbind(df2,meta[df2$Sample_ID,])

	write.csv(df2,file = paste0(prefix,"TC_sample_meta.csv"))

	##########################Exploring from here
	plotdata = scale(df2[,1:nsamples])
	png(paste0(prefix,"Sample_TC_all.png"),height = 600,width = 500)
		Heatmap(plotdata,row_split = df2[,featureOfInterest])
	dev.off()
	correlations = c()
	for(i in 1:nsamples){
		temp = cor.test(df2[,i],as.numeric(df2[,featureOfInterest]),method = cor_method)
		correlations = rbind(correlations,c(i,temp$estimate,temp$p.value))
	}
	colnames(correlations) = c("Comp","rho","p.value")
	correlations = data.frame(correlations)
	correlations$p.adj.fdr = p.adjust(correlations$p.value,method = "fdr")
	write.csv(correlations,file = paste0(prefix,"TC_sample_correlation.csv"))


	### TC2
	df_loading = data.frame(cp_decomp$U[[2]])
	colnames(df_loading) = paste0("Comp",1:ncol(df_loading))
	rownames(df_loading) = varFeatures
	write.csv(df_loading,file = paste0(prefix,"TC2_gene_loading.csv"))

	selectComp = correlations$p.adj.fdr<0.05
	if(sum(selectComp)>0){

		df_loading_Sig = df_loading[,correlations[selectComp,"Comp"]]
		rownames(df_loading_Sig) = varFeatures


		plotdata = scale(df_loading_Sig)
		select_genes = rowSums(abs(plotdata)>zthres)>0
		print(sum(select_genes))

		if(sum(select_genes)>0){
			png(paste0(prefix,"TC2_sig_gene_loadings.png"))
			print({
			Heatmap(plotdata[select_genes,])
			})
			dev.off()
		}
		select_genes2 = select_genes
		## TC1
		df_loading = data.frame(cp_decomp$U[[1]])
		colnames(df_loading) = paste0("Comp",1:ncol(df_loading))
		rownames(df_loading) = varFeatures
		write.csv(df_loading,file = paste0(prefix,"TC1_gene_loading.csv"))


		df_loading_Sig = df_loading[,correlations[selectComp,"Comp"]]
		rownames(df_loading_Sig) = varFeatures


		plotdata = scale(df_loading_Sig)
		select_genes = rowSums(abs(plotdata)>zthres)>0
		print(sum(select_genes))

		if(sum(select_genes)>0){
			png(paste0(prefix,"TC1_sig_gene_loadings.png"))
				print({
			Heatmap(plotdata[select_genes,])
			})
			dev.off()
		}
		select_genes1 = select_genes
		## Combined TC1, TC2
		select_genes = (select_genes2+select_genes1)>0
		write.csv(select_genes,file = paste0(prefix,"Selected_genes_from_TC1_TC2.csv"))
		################################################################################
		## Retrieve and combine negative and positive edges from the 
		## network for all samples within a category
		types = levels(meta[,featureOfCategory])
		nfeatures = length(varFeatures)
		net_pos = list()
		net_neg = list()
		net = list()
		for(i in types){
			net[[i]] = matrix(0,nfeatures,nfeatures)
			rownames(net[[i]]) = varFeatures
			colnames(net[[i]]) = varFeatures
			net_pos[[i]] = net[[i]]
			net_neg[[i]] = net[[i]]
			for(ssample in which(df2[,featureOfInterest] == i)){
				v1 = vec(t_coexp[,,ssample])
				m1 = matrix(v1,nrow  =  nfeatures)
				net[[i]] = net[[i]]+(m1>0)
				net_pos[[i]] = net_pos[[i]]+(m1>0)
				net_neg[[i]] = net_neg[[i]]-(m1<0)
			#m1 = m1[,colSums(m1)! = 0]
			}
		}

		## Collapse the negative and positive edges from each category by addition, 
		## and normalize by the total samples per category
		for(i in  types){
			net_sp = net_pos[[i]]+net_neg[[i]]
			for(j in 1:nrow(net_sp)){net_sp[j,j] = 0}
			net_sp = net_sp[select_genes, colSums( net_sp[select_genes,] ) >0 ]
			net_sp = net_sp/sum(df2[,featureOfInterest] == i)*100

			net_pair = c()
			for(j in 1:nrow(net_sp)){
				for( k in 1:ncol(net_sp)){
					if(net_sp[j,k]!= 0){
						net_pair = rbind(net_pair, c(rownames(net_sp)[j],colnames(net_sp)[k],net_sp[j,k]))
					}
				}
			
			}
			## save network for each category
			write.csv(net_pair,file = paste0(prefix,"net_",i,".csv"),quote = FALSE,row.names = F)
		}
		## save node annotation of the sig genes
		write.csv(cbind(select_genes[select_genes],"Sig"),file = paste0(prefix,"node.csv"))
	}
}

#' Extract the genes with highest loadings to TCs
#' @param cp_decomp decomposed matrices from HOSVD
#' @param selectComp selected TCs based on the correlation
#' @param varFeatures gene names
#' @param zthres threshold for z-score of loading values to identify significant driving genes per TC
#' @return A list of driving genes
selectGenesCo=function(cp_decomp,selectComp,varFeatures,zthres=2){
	df_loading=data.frame(cp_decomp$U[[1]])
	colnames(df_loading)=paste0("Comp",1:ncol(df_loading))
	rownames(df_loading)=varFeatures
	df_loading_Sig=df_loading[,selectComp,drop=FALSE]
	rownames(df_loading_Sig)=varFeatures
	plotdata=scale(df_loading_Sig)
	select_genes1=rowSums(abs(plotdata)>zthres)>0
	df_loading=data.frame(cp_decomp$U[[2]])
	colnames(df_loading)=paste0("Comp",1:ncol(df_loading))
	rownames(df_loading)=varFeatures
	df_loading_Sig=df_loading[,selectComp,drop=FALSE]
	rownames(df_loading_Sig)=varFeatures
	plotdata=scale(df_loading_Sig)
	select_genes2=rowSums(abs(plotdata)>zthres)>0
	select_genes = unique(c(names(select_genes2[select_genes2]),names(select_genes1[select_genes1])))
	return(select_genes)
}


#' Extract the tensor component (TC) that correlated with the feature of interest
#' @param cp_decomp decomposed matrices from HOSVD
#' @param nsamples number of samples
#' @param varFeatures gene names
#' @param meta metadata for samples
#' @param featureOfInterest column name of the feature of interest in the metadata: assuming the features are levels or numerical values
#' @param cor_method correlation method
#' @return A list of selected TCs
selectTC=function(meta,featureOfInterest,nsamples,cp_decomp,varFeatures,cor_method="spearman"){
	plotdata=scale(meta[,1:nsamples])
	correlations=c()
	for(i in 1:nsamples){
		temp=cor.test(meta[,i],as.numeric(meta[,featureOfInterest]),method=cor_method)
		correlations=rbind(correlations,c(i,temp$estimate,temp$p.value))
	}
	colnames(correlations)=c("Comp","rho","p.value")
	correlations=data.frame(correlations)
	correlations$p.adj.fdr=p.adjust(correlations$p.value,method="fdr")
	df_loading=data.frame(cp_decomp$U[[2]])
	colnames(df_loading)=paste0("Comp",1:ncol(df_loading))
	rownames(df_loading)=varFeatures
	selectComp=correlations$p.adj.fdr<0.05
	return(selectComp)
}


#' Permutation on samples to extract the genes
#' @param nsamples number of samples
#' @param varFeatures gene names
#' @param meta metadata for samples
#' @param featureOfInterest column name of the feature of interest in the metadata assuming the features are levels or numerical values
#' @param cor_method correlation method
#' @param zthres threshold for z-score of loading values to identify significant driving genes per TC
#' @export
permTest = function(t_coexp, meta,varFeatures,  featureOfInterest, featureOfCategory, prefix, nsamples_per_group = 5){
	if(!dir.exists(file.path(prefix))){
		prefix = "."
	}
	dir.create(file.path(prefix, "perm"), showWarnings  =  FALSE)
	prefix = paste0(prefix,"/perm/")
	conditions = levels(featureOfCategory)
	sig_Genes=foreach(j=1:100,.combine=rbind)%dopar% {
		sample_sub=c()
		for(i in conditions){		
			sample_sub=c(sample_sub,sample(which(meta[,featureOfInterest]==i),nsamples_per_group,replace=sum(meta[,featureOfInterest]==i)<nsamples_per_group))
		}
		meta_sub=meta[sample_sub,]
		nsamples=length(sample_sub)
		
		cp_decomp <- rTensor::cp(t_coexp[,,sample_sub], num_components = nsamples, max_iter = 100)

		df2=data.frame(cp_decomp$U[[3]])
		colnames(df2)=paste0("Comp",1:ncol(df2))
		df2$Sample_ID=sampleNames[sample_sub]
		df2=cbind(df2,meta_sub[df2$Sample_ID,])

		selectComp=selectTC(df2,nsamples,cp_decomp,varFeatures)
		if(sum(selectComp)>0){
			select_genes=selectGenesCo(cp_decomp,selectComp,varFeatures)
			return(cbind(select_genes,j))
		}else{return(NULL)}
	}
	write.csv(sig_Genes,file=paste0(prefix,"perm_",nsamples_per_group,".csv"))
}