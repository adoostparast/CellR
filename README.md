# CellR
CellR is a single cell-based data-driven method to recover and quantify the cellular composition of bulk transcriptional data. It is a fully unsupervised approach based on clustering the reference single-cell RNA-Seq (scRNA-seq) followed by extracting the unique marker genes defining each cell cluster. 
# Arguments
CellR receives seven input arguments as follows:
* **Reference scRNA-seq data**
Reference scRNA-seq data is a tab delimited text file whose rows represent the unique gene names and each column represents a cell. Note that the upper left most element of the fule should be empty.
CellR receives scRNA-seq file as non-normalized raw counts.
* **Bulk RNA-seq data**
