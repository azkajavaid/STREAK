# STREAK 
<img width="1120" alt="Screenshot 2023-03-14 at 8 40 02 AM" src="https://user-images.githubusercontent.com/46649424/225003703-683fedcb-115a-47e3-aefb-12264169fb18.png">

### Abstract
The accurate estimation of cell surface receptor abundance for single cell transcriptomics data is important for the tasks of cell type and phenotype categorization and cell-cell interaction quantification. We previously developed an unsupervised receptor abundance estimation technique named SPECK (Surface Protein abundance Estimation using CKmeans-based clustered thresholding) to address the challenges associated with accurate abundance estimation. In that paper, we concluded that SPECK results in improved concordance with Cellular Indexing of Transcriptomes and Epitopes by Sequencing (CITE-seq) data relative to comparative unsupervised abundance estimation techniques using only single-cell RNA-sequencing (scRNA-seq) data. In this paper, we outline a new supervised receptor abundance estimation method called STREAK (gene Set Testing-based Receptor abundance Estimation using Adjusted distances and cKmeans thresholding) that leverages associations learned from joint scRNA-seq/CITE-seq training data and a thresholded gene set scoring mechanism to estimate receptor abundance for scRNA-seq target data. We evaluate STREAK relative to both unsupervised and supervised receptor abundance estimation techniques using two evaluation approaches on three joint scRNA-seq/CITE-seq datasets that represent two human tissue types. We conclude that STREAK outperforms other abundance estimation strategies and provides a more biologically interpretable and transparent statistical model.
