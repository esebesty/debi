
============================================================

README file for DeBi algorithm 
============================================================

Reference:

Akdes Serin and Martin Vingron DeBi: Discovering Differentially Expressed
Biclusters using a Frequent Itemset Approach Algorithms Mol Biol. 2011; 6: 18.
doi:  10.1186/1748-7188-6-18

============================================================

Installation
============================
cd src
make


Parameters
============================
The name of the executable is debi and the parameters are set as 

./debi infile outFilePath

inFile: is the tab deliminated file with the first row 'sample ids' and first column 'gene ids'. (see the samples folder for example inputs)  outFilePath: is the path of the biclustering results folder

-b binarization level (for binary data not needed otherwise it must be defined) 
-o maximum overlap level of the biclusters default=0.5 (overlap=1: no overlap is allowed, overlap=0: all the biclusters are taken, no filtering applied)
-p pattern of regulation up or down (u, d) default=u

output file format
<bicluster gene size> <bicluster sample size> <bicluster score> <normalized bicluster score>
genes in bicluster 
samples in bicluster


Sample Commands
============================

---yeast data---

up regulated biclusters (binarization level=log10(1.58))

./debi "../samples/yeast/hugues.txt" "../results/" -b0.1986571 -pu 

down regulated biclusters (binarization level=-log10(1.58))

./debi "../samples/yeast/hugues.txt" "../results/" -b-0.1986571 -pd

---MSIGDB data--- 

first download the data from http://www.broadinstitute.org/gsea/msigdb/ make the binary matrix and put it into samples/msigdb folder 

./debi "../samples/msigdb/Hs.gmtl.c2.txt" "../results/" -b1 -pu


---DLBLC data---

./debi "../samples/dlbcl/exp_dlbcl.txt" "../results/" -b1 -pu
./debi "../samples/dlbcl/exp_dlbcl.txt" "../results/" -b-1 -pd


---ExPO data--- 

first download the data from  'https://expo.intgen.org/geo/home.do log2 transform, quantile normalize and mean subtract and put it into samples/cancer_data folder 

./debi "../samples/expo/expO_quantile_log2_mean_WH.txt" "../results/" -b1  -pu 
./debi "../samples/expo/expO_quantile_log2_mean_WH.txt" "../results/"  -b-1  -pd


---CMAP data---

first download the 'ftp://ftp.broad.mit.edu/pub/cmap/ratioMatrix.txt' data and put it into samples/cmap folder 

./debi "../samples/cmap/ratioMatrix.txt" "../results/"  -b2  -pu 
./debi "../samples/cmap/ratioMatrix.txt" "../results/" -b0.5  -pd 


