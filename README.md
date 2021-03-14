# Non-linear-and-linear-techniques-for-dimensionality-reduction-visualization-on-single-cell-data.
This implementation is part of my thesis.

## Contents
1. Introduction
2. Getting started
3. Datasets Used
4. Installation
5. Running the code
6. Contributing to the project
7. Acknowledgments
8. License

## 1. Introduction
Single-cell sequencing (scRNA-seq) is an emerging technology used to capture cell information at a single-nucleotide resolution and by which individual cells can be analyzed separately. Single-cell data is high-dimensional and sparse data lead to some analytical challenges. Analyzing scRNA-seq data can be divided into two main categories: at the cell level and gene level. Finding cell sub-networks or highly deferentially expressed tissue-specific gene lists is one of the common challenges at the cell level.

Grouping cells into different clusters to find heterogeneity is one of the significant steps in single-cell data analysis. However, due to high-dimensional data its uncertain to get good clustering and visualization. Hence, non-linear dimensionality reduction techniques such as MLLE are efficient, and linear methods like ICA are excellent in visualization. The combination of both techniques combined with clustering gives the best clustering scores. 

## 2. Getting Started
The whole code is written in Python (3.4+).
Most parts of the implementation are using Scanpy package, used for single-cell analysis.
You need below packages to replicate this work.
- [Scanpy](https://scanpy.readthedocs.io/en/stable/)
- Numpy
- Pandas
- [Scikit-learn](https://scikit-learn.org/stable/) – metrics, KMeans, FastICA, LocallyLinearEmbedding, Isomap
- [Plotly.graph_objects](https://plotly.com/python/graph-objects/) 

## 3. Datasets used
Single-cell RNA sequencing data.
| Dataset Name  | Accession Number | Number of Datasets |
| ------------- | ------------- | ------------- |
| SARS-Cov  | [GSE148729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148729)  | 2 |
| Muraro  | [GSE85241](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241_cellsystems_dataset_4donors_updated.csv.gz)  | 1 |
| Segerstolpe  | [E_MTAB_5061](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt)  | 1 |
| Xin  | [GSE81608](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81608/suppl/GSE81608_human_islets_rpkm.txt.gz)  | 1 |
| Wang  | [GSE83139](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83139/suppl/GSE83139_tbx-v-f-norm-ntv-cpms.csv.gz)  | 1 |
| PBMC  | [10X Genomics](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)  | 1 |
| Baron  | [GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133)  | 6 |

## 4. Installation
1.	Clone the repository.
```
git clone https://github.com/saiteja-danda/Non-linear-and-linear-techniques-for-dimensionality-reduction-visualization-on-single-cell-data.git
```
3.	Install necessary libraries as mentioned.
4.	Download the data using accession numbers.
5.	Install Spyder or Jupyter notebook or any other python IDE to run the experiment.

## 5. Running the code
Input can be read in the form of .txt, .csv and .mtx files. All the datasets mentioned above are in these formats only.
You can easily implement this project as every cell in the code file has comments and describes why that particular piece of code is used.

## 6. Contributing to the project
Open source communities are such unique places to learn, innovate, and contribute. Any contributions to this project are deeply appreciated.

## 7. Acknowledgments
This research was partially supported by Mitacs and Natural Sciences and Engineering Research Council of Canada, NSERC. I want to thank [Dr. Luis Rueda](https://luisrueda.myweb.cs.uwindsor.ca/) for his continuous support and motivation, Akram Vasighizaker, a PhD student for her collaboration on this project, and the University of Windsor Office of Research and Innovation.

## 8. License
  See [License](https://github.com/saiteja-danda/Non-linear-and-linear-techniques-for-dimensionality-reduction-visualization-on-single-cell-data./blob/main/LICENSE)
