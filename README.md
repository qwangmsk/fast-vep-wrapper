fast-vep-wrapper
================

This project aims to use VEP to quickly and conveniently annotate variants stored in large MAF files. Typical scenario to run fast-vep-wrapper.pl is to annotate MAF files containing somatic variants from thousands of or more tumor samples. 

Installation
------------

The script requires users to install VEP and [vcf2maf](https://github.com/ckandoth/vcf2maf).

To improve portability, it utilizes a configuration file to store paths of vep, ref_fasta, and vcf2maf, together with other parameters. By default, the script uses a file named config.txt in the same directory as fast-vep-wrapper.pl as configuration file, unless user specifies a configuration file explicitly in command line.  

Description
-----------

[fast-vep-wrapper.pl]() extracts unique variants from input MAF file. Here, 'unique' means unique combination of Chromosome, Start_Position, Reference_Allele and Tumor_Allele. If an annotated MAF file is provided, variants in the annotated MAF will not be re-annotated. Only new variants in the MAF file are processed (using [vcf2maf](https://github.com/ckandoth/vcf2maf)). For new variants, 'TUMOR' and 'NORMAL' are used as'Tumor_Sample_Barcode' and 'Matched_Norm_Sample_Barcode' to reduce multithreading overhead.


Another script file [run-vep-wrapper.pl]() provides application-level wrapper functions. I mainly use it to test [fast-vep-wrapper.pl](). 

Acknowledgements
----------------

Thank Cyriac Kandoth, Frederick Criscuolo, Onur Sumer, and Sumit Middha for insightful suggestions and helpful discussion

Authors
-------

Qingguo Wang (josephw10000@gmail.com)

