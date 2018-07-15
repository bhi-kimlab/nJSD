#nJSD

nJSD is code for calculating distance between two biological networks instantiated with gene-expression profiles using entropy concept. 
It was designed to measure intratumor heterogeneity from bulk RNA-sequencing data.
Transcriptome-based ITH (tITH) of tumor state was calculated by considering both normal state and ideally heterogeneous state.

Requirements
---------------------
* Linux
* Python 2.7+    (recommend 2.7.15)
* NumPy 
* NetworkX 2.1+  (recommend 2.1)


 nJSD
----------------------
It compute distance (nJSD) between two GEPs and calculate tITH score.
Two functions are supported.
 'whole' is to caclulate nJSD on whole given gene set in expression data.
 'geneset' need additional geneset list file and show tITH on user given gene set.
   
    nJSD.py whole [-h] [-n NETWORK] [-r R_GEP] [-i Q_GEP]
    nJSD.py geneset [-h] [-n NETWORK] [-r R_GEP] [-i Q_GEP] [-t GENESET_FILE]
    
    -h, --help  show this help message and exit
    -n NETWORK  Location to network file: geneA geneB
    -r R_GEP    File name of Refernece gene-expression profile
    -i Q_GEP    File name of Query gene-expression profile
    -t GENESET  File name of Gene set (like pathway)

Network file must follow below format.

    GeneA GeneB               # Header
    GeneSymbol1 GeneSymbol2
    GeneSymbol1 GeneSymbol3
    GeneSymbol1 GeneSymbol4
    ...

GEP file must follow below format.

    GeneSymbol  ExpressionValue       # Header
    GeneA 10
    GeneB 20
    BeneC 30
    ...

Geneset List file must follow below format. (No header)

    Group1Name  GeneA   GeneB   GeneC   ...
    Group2Name  GeneD   GeneE   GeneF   ...
    Group3Name  GeneA   GeneG   GeneH   ...
    ...


When the gene set of reference GEP is differ to gene set of query GEP file and geneset list file.
The difference is dumped into a file with name "dumpgene+date".


Toy Data
----------------------

In the example directory, there are toy data.

example:

    python nJSD.py whole -n example/Toy.network -r example/Toy.profile1 -i example/Toy.profile2
    
result:
    
    example/Toy.profile2    [Ref. -> Query: 0.003935]       [Query -> stateH: 0.007693]     <tITH: 0.338413>
    
example:

    python nJSD.py geneset -n example/Toy.network -r example/Toy.profile1 -i example/Toy.profile2 -t example/Toy.geneset

result:

    example/Toy.profile2     1st_pwy         [Ref. -> Query: 0.007822]      [Query -> stateH: 0.009385]     <tITH: 0.454582> 
    example/Toy.profile2     3rd_pwy         [Ref. -> Query: 0.005215]      [Query -> stateH: 0.007102]     <tITH: 0.423379> 
    example/Toy.profile2     2nd_pwy         [Ref. -> Query: 0.000000]      [Query -> stateH: 0.004261]     <tITH: 0.000000> 
    example/Toy.profile2     4th_pwy         [Ref. -> Query: 0.007909]      [Query -> stateH: 0.007403]     <tITH: 0.516523> 
    example/Toy.profile2     5th_pwy         [Ref. -> Query: 0.004470]      [Query -> stateH: 0.012175]     <tITH: 0.268536> 


Citation
----------------------
Y. Park, S. Lim, J. Nam, S. Kim, Measuring intratumor heterogeneity by network entropy using RNA-seq data, Scientific Reports (2016)
