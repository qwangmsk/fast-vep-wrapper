fast-vep-wrapper
================

This project aims to use VEP to quickly and conveniently annotate mutations stored in large MAF files. Typical scenario to run [fast-vep-wrapper.pl]() is to annotate MAF files containing somatic mutations from thousands of or more tumor samples. 

Installation
------------

The script requires users to install [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) and [vcf2maf](https://github.com/ckandoth/vcf2maf).

To improve portability, it supports (but does not require) the use of a configuration file to store paths of vep, ref_fasta, and [vcf2maf](https://github.com/ckandoth/vcf2maf), together with other parameters. In the configuration file, parameters are defined as "variable=value" pairs. By default, the script uses a file named config.txt in the same directory as [fast-vep-wrapper.pl]() as configuration file, unless user specifies a different file explicitly in command line.

Description
-----------

[fast-vep-wrapper.pl]() extracts unique variants from input MAF file. Here, 'unique' means unique combination of Chromosome, Start_Position, Reference_Allele and Tumor_Allele. If an annotated MAF file is provided, variants in the annotated MAF will not be re-annotated. Only new variants in the MAF file are processed (using [vcf2maf](https://github.com/ckandoth/vcf2maf)). New variants are annotated in a way to minimize multithreading overhead in order to improve speed.

[fast-vep-wrapper.pl]() is derived from [vcf2maf](https://github.com/ckandoth/vcf2maf). By re-using previous annotation and improving parallel efficiency, it enhances annotation efficiency that is lacking in [vcf2maf](https://github.com/ckandoth/vcf2maf). The speed of [fast-vep-wrapper.pl]() is linear to the number of CPUs (vep_forks). For large MAF files and when using multiple CPUs, it can easily reduce computation time from hours to minutes.

Besides [fast-vep-wrapper.pl](), this project also consists of the following files:
 1. Script file [run-vep-wrapper.pl]() provides application-level wrapper functions. I mainly use it to test [fast-vep-wrapper.pl](). I also use [run-vep-wrapper.pl]() to run [maf2maf.pl](https://github.com/ckandoth/vcf2maf). Unlike [fast-vep-wrapper.pl](), [run-vep-wrapper.pl]() requires users to define a configuration file.
 2. Text file [depth_cols.txt]() contains read count-related column names (1st column) that appear repeatitively in various MAF files. The 2nd column of the file are [vcf2maf](https://github.com/ckandoth/vcf2maf) arguments corresponding to these count-related column names. When parsing MAF file header, [run-vep-wrapper.pl]() looks for words in the 1st column of [depth_cols.txt]() and then automatically constructs correct command to run [fast-vep-wrapper.pl]() / [maf2maf.pl](https://github.com/ckandoth/vcf2maf). The file [depth_cols.txt]() allows fully automatic annotation of MAF files from diverse sources.
 3. Script [msk-vep-wrapper.pl]() stores and re-uses VEP-annotated MAF files. It is designed specifically for MSKCC user who does not want to manage previous annotation files.

Acknowledgements
----------------

    Thank Cyriac Kandoth, Frederick Criscuolo, Onur Sumer, Sumit Middha, and Benjamin Gross for helpful suggestions and discussion

Authors
-------

    Qingguo Wang
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065

