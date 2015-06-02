fast-vep-wrapper
================

This project aims to use VEP to quickly annotate variants stored in large MAF files. Typical scenario to run fast-vep-wrapper.pl is to annotate MAF files containing somatic variants from thousands of or more tumor samples. 

Installation
------------

The script requires users to install VEP and vcf2maf (https://github.com/ckandoth/vcf2maf).

To improve portability, it utilizes a configuration file to store paths of vep, ref_fasta, and vcf2maf, together with other parameters. In default settings (config.txt as default file name; the same path as this program as default directory to store the configuration file), user does not need to specify the configuration file in command line.  

Method
------

The script extracts unique variants from input MAF file. Here, 'unique' means unique combination of Chromosome, Start_Position, Reference_Allele and Tumor_Allele. If an annotated MAF file is provided, variants in the annotated MAF will not be re-annotated. Only new variants in the MAF file are processed (using vcf2maf: https://github.com/ckandoth/vcf2maf). For new variants, 'TUMOR' and 'NORMAL' are used as'Tumor_Sample_Barcode' and 'Matched_Norm_Sample_Barcode' to reduce multithreading overhead.


Derived from vcf2maf, this program is designed to enhance annotation efficiency that is lacking in vcf2maf. By reusing previous annotation and improving parallel efficiency, it reduces time for annotating large MAF files from hours to minutes.

Acknowledgements
----------------

Thank Sumit Middha, Frederick Criscuolo, and Cyriac Kandoth for suggestion and discussion

Authors
-------

Qingguo Wang (josephw10000@gmail.com)

