# CoSig
This R package is designed for generating co-expression network for single cell dataset and extracting gene signatures that reflect network variations along features of interest.

# Installation
devtools::install_github("bingqing-Xie/CoSig")

# Usage
```
## Proceed with any seurat pipeline to obtain a seurat object (obj).
## Set parallel threads
initializeMultiCores(cores=10)
## Network will be saved in the prefix/raw/ folder
results = computeCoExp(obj, prefix = ".", sampleid = "Sample_ID",
                        recompute_glasso = TRUE, recompute_HOSVD = TRUE, ReComputeVarGene = FALSE,
                        mean_thres = 0.1)
cp_decomp = results[[1]]
t_coexp = results[[2]]
sampleNames = results[[3]]
varFeatures = results[[4]]

## Analysis results will be saved in the prefix/result/ folder
meta = read.csv("PATH/TO/YOUR/META.csv")
featureOfInterest = "YOUR_FEATURE_IN_META"  # Numerical features to be analyzed
featureOfCategory = "YOUR_FEATURE_IN_META"  # This can be the same feature as above, and needs to be a factor
featureAnalysis(prefix = "." , cp_decomp, t_coexp,
                sampleNames,varFeatures, meta,
                featureOfInterest, featureOfCategory,
                cor_method  =  "spearman", zthres = 2)

## If permutation test is needed, results will be saved in the prefix/perm/ folder
## for randomly selected samples per featureOfCategory group
permTest(t_coexp, meta,varFeatures,  featureOfInterest, featureOfCategory, prefix, nsamples_per_group = 5)
```
