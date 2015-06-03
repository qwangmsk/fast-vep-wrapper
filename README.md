fast-vep-wrapper
================

This project aims to use VEP to quickly and conveniently annotate mutations stored in large MAF files. Typical scenario to run [fast-vep-wrapper.pl]() is to annotate MAF files containing somatic mutations from thousands of or more tumor samples. 

Installation
------------

The script requires users to install [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) and [vcf2maf](https://github.com/ckandoth/vcf2maf).

To improve portability, it utilizes a configuration file to store paths of vep, ref_fasta, and [vcf2maf](https://github.com/ckandoth/vcf2maf), together with other parameters. In configuration file, parameters are defined as "variable=value" pairs. By default, the script uses a file named config.txt in the same directory as [fast-vep-wrapper.pl]() as configuration file, unless user specifies a configuration file explicitly in command line.  

Description
-----------

[fast-vep-wrapper.pl]() extracts unique variants from input MAF file. Here, 'unique' means unique combination of Chromosome, Start_Position, Reference_Allele and Tumor_Allele. If an annotated MAF file is provided, variants in the annotated MAF will not be re-annotated. Only new variants in the MAF file are processed (using [vcf2maf](https://github.com/ckandoth/vcf2maf)). For new variants, 'TUMOR' and 'NORMAL' are used as'Tumor_Sample_Barcode' and 'Matched_Norm_Sample_Barcode' to reduce multithreading overhead.

[fast-vep-wrapper.pl]() is derived from [vcf2maf](https://github.com/ckandoth/vcf2maf). By re-using previous annotation and improving parallel efficiency, it enhances annotation efficiency that is lacking in [vcf2maf](https://github.com/ckandoth/vcf2maf). The speed of [fast-vep-wrapper.pl]() is linear to the number of CPUs (or vep_forks). For large MAF files and when using multiple CPUs, it can easily reduce computation time from hours to minutes.

Another script file [run-vep-wrapper.pl]() provides application-level wrapper functions. I mainly use it to test [fast-vep-wrapper.pl](). 

Acknowledgements
----------------

Thank Cyriac Kandoth, Frederick Criscuolo, Onur Sumer, Sumit Middha, and Benjamin Gross for helpful suggestions and discussion

Authors
-------

Qingguo Wang

Nikolaus Schultz Lab

Memorial Sloan Kettering Cancer Center

New York, NY 10065

