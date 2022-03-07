# biopython_gbk_parser
A GBK parser using biopython, which will provide all loci in a genome given a GBK file of a genomic assembly. The second argument is the name of a locus tag in said genome, which the parser will use to print the genomic neighbourhood of said locus tag.

How to use:\
biopython_gbk_parser.py  input_file  locus

### Arguments:
- input_file: the GBK file, (an example of a .gbff file is available in this repository)
- locus: the name of a locus tag in the GBK file, it's prior knowledge of the locus of interest.
  
 
### Functions:
  - make_list_and_dict(input_file): it will parse the GBK file to provide a list_of_genes (which contains all loci in order) and a dictionary dict_of_genes (locus : [strand, biological product])
  - get_genomic_neighbourhood_list(locus_tag,list_of_genes): it will provide a genomic_neighbourhood_list of the loci -2,-1,+1 y +2 surrounding our locus_tag we use as query, including the locus_tag itself
  - print_result(genomic_neighbourhood_list,dict_of_genes): it will print the genomic neighbourhood using the format method, each locus will be printed in one direction or another depending on the strand the locus is located at:
    - positive strand: [ locus ]>
    - negative strand: <[ locus ]
