#nJSD

nJSD is code for calculating distance between two biological networks instantiated with gene-expression profiles using entropy concept. 
It was designed to measure intratumor heterogeneity from bulk RNA-sequencing data.
Transcriptome-based ITH (tITH) of tumor state was calculated by considering both normal state and ideally heterogeneous state.

Requirements
---------------------
* Linux/Unix
* Python 2.7
* NumPy 
* NetworkX


Simple nJSD
----------------------
It compute distance (nJSD) between two GEPs and calculate tITH score.


    run_simple.py [-h] [-n NETWORK] [-r R_GEP] [-q Q_GEP]
    -h, --help  show this help message and exit
    -n NETWORK  Location to network file: geneA geneB
    -r R_GEP    File name of Refernece gene-expression profile
    -q Q_GEP    File name of Query gene-expression profile


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
    ....



In the example directory, there are test data.

example:

    ./run_simple.py -n example/Toy.network -r example/Toy.profile1 -q example/Toy.profile2


result:

    example/Toy.profile2    [Ref. -> Query: 0.003935]       [Query -> stateH: 0.007693]     <tITH: 0.338413>




Citation
----------------------
Y. Park, S. Lim, J. Nam, S. Kim, Measuring intratumor heterogeneity by network entropy using RNA-seq data, Scientific Reports (in press)
