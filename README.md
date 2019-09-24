![Cellr](https://user-images.githubusercontent.com/24727526/65530825-aaaf0b00-dec6-11e9-9531-ca8a72a9abd4.JPG)

# CellR
CellR is a single cell-based data-driven method to recover and quantify the cellular composition of bulk transcriptional data. It is a fully unsupervised approach based on clustering the reference single-cell RNA-Seq (scRNA-seq) followed by extracting the unique marker genes defining each cell cluster. 


# Installation
CellR bears some dependencies including dplyr, Matrix, Seurat (v 2.3.0, To install, check this out: https://satijalab.org/seurat/install.html), edgeR, glmnet. Please make sure to have devtools installed on your machine. CellR automatically installs the required dependencies during the initial installation. However, in case of failure, users can manually install the mentioned packages. CellR is built under R vesrion 3.5.1.

CellR can be installed from source after downloading the CellR zip file using: install.packages("CellR_0.1.0.tar.gz"). Next, it can be loaded using: library('CellR').
# Arguments
The general format for calling CellR is: 
Output<-Deconvolution(RNA-seq, scRNA-seq, GTEx, Cell-types, Minimum-cells, Minimum-genes, Dimension, Alpha, Maximum_gene_counts, Minimum_gene_counts)

CellR arguments are explined in detail as follows:
* **Reference scRNA-seq data**
Reference scRNA-seq data is a tab delimited text file whose rows represent the unique gene symbols and each column represents a cell. Note that the upper left most element of the file should be empty.
CellR receives scRNA-seq file as a non-normalized raw counts matrix.

* **Bulk RNA-seq data**
Bulk RNA-seq data is a tab delimited text file whose rows represent the unique gene symbols and each column represents a sample. Note that the upper left most element of the file should be empty.
CellR receives RNA-seq file as a non-normalized raw counts matrix.

* **GTEx TPM data**
GTEx data is a tab delimited text file whose rows represent the unique gene symbols and each column represents an individual sample. Note that the upper left most element of the file should be empty. For any human tissue, relevant GTEx data can be obtained from GTEx portal: https://gtexportal.org/home/.

* **Cell-types**
A column vector with rows indicating the cell IDs exactly as is in the reference scRNA-seq data, and the column denotes the cell type of each cell. Note that the first cell in the first column should be empty.

* **Minimum-cells**
CellR keeps all the genes being expressed in > Minimum-cells (Minimum-cells>=3 is used in our experiments).

* **Minimum-genes**
CellR keeps all the cells with at least "Minimum-genes" identified genes (Minimum-genes>=200 is used in our experiments).

* **Dimension**
The total number of the top pincipal components (PCs) used for clustering the reference scRNA-seq data. Clustering is the backbone of CellR, so we strongly recommend accurately setting Dimension. We used Dimension=12 in our manuscript. Readers are encouraged to read more about how to obtain the best number of PCs for clustering: https://satijalab.org/seurat/pbmc3k_tutorial.html

* **Alpha**
Alpha is a term used for taking into account inter/intra cellular differences during the optimization process. Alpha=1 indicates lasso mode and Alpha=0 denotes the ridge mode.

* **Maximum_gene_counts**
Genes with unique counts over 'Maximum_gene_counts' will be filtered out.

* **Minimum_gene_counts**
Genes with unique counts less than 'Minimum_gene_counts' will be filtered out.

# A real example to run
CellR bears some 
https://www.dropbox.com/sh/pikqwjdf0bs2t88/AAAlHVe4B0pKagtTU_3BuFvva?dl=0

# Output
CellR outputs the percentage of the proportion of each identified cluster from the reference scRNA-seq data within the bulk RNA-seq sample as well as the identified cell type-specific markers used during the deconvolution process.
