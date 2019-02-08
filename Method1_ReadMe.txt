This method finds all core genes that contain 10-40 (number of progenies, by default) values below a cut-off in at least LOW_BOUND consecutive windows. Less genes are detected by this method, compared to Method2.

1. To find deletion windows, type:
python3 find_deletion_windows.py read.data.txt 

The output is 'deletion_windows.txt'.

There are three parameters you can change:
-c: cut_off, the default is 0.01
-l: low_bound, which determines the number of consecutive windows. The default is 3 
-i: interval, which restrain the number of deletions in a picked window. 
The default is [10, 40]. Takes two integers. Example input: 10 40


2. To find deleted genes in those windows, type:
python3 find_deleted_genes.py deletion_windows.txt PlasmoDB-39_Pfalciparum3D7.gff 

The output is 'deleted_genes.txt'.


3. To find deleted genes in the core genome, type:
python3 find_deleted_core_genes.py deleted_genes.txt core.txt 

The output is 'deleted_core_genes.txt'.
