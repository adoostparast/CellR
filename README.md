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
For convenient use of CellR, we have created a data repository which contains the data on Alzheimer's bulk RNA-seq data being used in the paper. In the following link, https://www.dropbox.com/sh/pikqwjdf0bs2t88/AAAlHVe4B0pKagtTU_3BuFvva?dl=0, you will find four data files including: Bulk.txt which is the RNA-seq data to be deconvolved, FCortex: the reference single-cell RNA-seq data from prefrontal cortex which is used by CellR to recover the marker genes, GTExExpessionSymbols.txt: the GTEx RNA-seq data used in CellR, Cells.txt: the identity of the cells in the FCortex.txt data. We encourage the users to check the structure of each data for smooth use of the toolkit in the future.
To run CellR using the supplied data, please use the following commands:

library(CellR)

Bulk<-read.table("Bulk.txt",header=TRUE)

Single<-read.table("FCortex.txt",header=TRUE)

GTEx<-read.table("GTExExpressionSymbol.txt",header=TRUE)

Cells<-read.table("Cells.txt",header=TRUE)

Output<-Deconvolution(Bulk,Single,GTEx,Cells,3,200,12,1,2500,1)

'Output' is an R object with three categories. Users can extract the results as follows:

Proportions<-Output$Proportion      # 'Proportions' yields a table whose rows are cell-types and the columns represent a bulk samples. Each cell within the table denotes the percentage of its corresponding cell-type within its corresponding sample. Below is a snapshot of this table.

![Capture](https://user-images.githubusercontent.com/24727526/65543165-a1ca3380-dede-11e9-82ab-b3ae398a3e90.JPG)


Markers<-Output$Markers   # 'Markers' is a table which shows the identified markers of each cell population which then are used in CellR. It contains the p-values of the difference of the gene expression of each signature gene in its corresponding cell cluster compared to others. Below is a snapshot of this table.

![Capture](https://user-images.githubusercontent.com/24727526/65543507-464c7580-dedf-11e9-9a42-e34b34910db5.JPG)

Clusters<-Output$Clusters  # 'Clusters' is a table showing the type (cluster) of each cell within the reference scRNA-seq data.




Note: Use of these arguments depends on the charachteristics of the data and the users can change these numbers in other applications.

------------------------------------------------------------------------------------------------------------------------
If you have any questions, please contact me at doostparaa@email.chop.edu.



# Output
CellR outputs an R object which includes the percentage of the proportion of each identified cluster from the reference scRNA-seq data within the bulk RNA-seq sample as well as the identified cell type-specific markers used during the deconvolution process.
