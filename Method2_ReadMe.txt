Method2 adopted methods used by Connor Howington in VCF.
This method finds all core genes that show at least LOW_BOUND zeros in at least one window. More genes are detected by this method, compared to Method1.

1. To find deletion windows, type:
python3 CNV_Match.py read.data.txt 

The output is 'data.map'.

There are three parameters you can change:
-l: low_bound, only return windows that have at least LOW_BOUND zeros. The default is 5 


2. To find deletion windows that are in core genome, type:
python3 CNV_in_core.py core.txt data.map 

The output is 'data.map.core'.


3. To find deleted core genes, type:
python3 CNV_in_gene.py data.map.core PlasmoDB-39_Pfalciparum3D7.gff 

The output is 'deleted_core_genes.txt'.
