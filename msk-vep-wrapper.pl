#!/usr/bin/env perl

# msk-vep-wrapper.pl - This script stores and re-uses VEP-annotated MAF files


use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Cwd 'abs_path';

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

my ( $maf2maf, $vep_wrap_script ) = ( '/ssd-data/cmo/opt/vcf2maf/maf2maf.pl', '' );
my ( $input_maf, $output_maf, $tmp_dir, $custom_enst_file );
my ( $tum_depth_col, $tum_rad_col, $tum_vad_col );
my ( $nrm_depth_col, $nrm_rad_col, $nrm_vad_col );
my ( $retain_cols, $vep_path, $vep_data, $vep_forks, $ref_fasta );
my ( $man, $help ) = ( 0, 0 );

GetOptions(
'help!' => \$help,
'man!' => \$man,
'input-maf=s' => \$input_maf,
'output-maf=s' => \$output_maf,
'tmp-dir=s' => \$tmp_dir,
'tum-depth-col=s' => \$tum_depth_col,
'tum-rad-col=s' => \$tum_rad_col,
'tum-vad-col=s' => \$tum_vad_col,
'nrm-depth-col=s' => \$nrm_depth_col,
'nrm-rad-col=s' => \$nrm_rad_col,
'nrm-vad-col=s' => \$nrm_vad_col,
'retain-cols=s' => \$retain_cols,
'custom-enst=s' => \$custom_enst_file,
'vep-path=s' => \$vep_path,
'vep-data=s' => \$vep_data,
'vep-forks=s' => \$vep_forks,
'ref-fasta=s' => \$ref_fasta,
'maf2maf=s' => \$maf2maf,
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );


# Check command line input
( defined $input_maf ) or die "ERROR: input MAF is not provided!\n";
( defined $maf2maf && -e $maf2maf ) or die "ERROR: maf2maf.pl is not specified\n";
( defined $tmp_dir && -d $tmp_dir) or die "ERROR: temporary dir is not specified or is not a directory\n";

my ( $maf_name, $maf_path ) = fileparse( abs_path( $input_maf ) );
# Construct path of cached VEP annotation file
my $annotated_maf = $maf_path;
$annotated_maf =~ s/\//_/g;
$annotated_maf .= 'TUMOR_vs_NORMAL.vep.maf';
$annotated_maf  = "$tmp_dir/$annotated_maf";

# Look for fast-vep-wrapper.pl in the same directory as this script
$vep_wrap_script = "$FindBin::Bin/fast-vep-wrapper.pl";
( -s $vep_wrap_script ) or die "Error: $vep_wrap_script does not exist\n";

# Create command string to run fast-vep-wrapper.pl
my $wrapper_cmd = "perl $vep_wrap_script --maf2maf $maf2maf --input-maf $input_maf --output-maf $output_maf";
$wrapper_cmd .= " --tum-depth-col $tum_depth_col"   if ( $tum_depth_col );
$wrapper_cmd .= " --tum-rad-col $tum_rad_col"       if ( $tum_rad_col );
$wrapper_cmd .= " --tum-vad-col $tum_vad_col"       if ( $tum_vad_col );
$wrapper_cmd .= " --nrm-depth-col $nrm_depth_col"   if ( $nrm_depth_col );
$wrapper_cmd .= " --nrm-rad-col $nrm_rad_col"       if ( $nrm_rad_col );
$wrapper_cmd .= " --nrm-vad-col $nrm_vad_col"       if ( $nrm_vad_col );
$wrapper_cmd .= " --retain-cols \"$retain_cols\""   if ( $retain_cols );
$wrapper_cmd .= " --custom-enst $custom_enst_file"  if ( $custom_enst_file );
$wrapper_cmd .= " --vep-path $vep_path"             if ( $vep_path );
$wrapper_cmd .= " --vep-data $vep_data"             if ( $vep_data );
$wrapper_cmd .= " --vep-forks $vep_forks"           if ( $vep_forks );
$wrapper_cmd .= " --ref-fasta $ref_fasta"           if ( $ref_fasta );
#$wrapper_cmd .= " --annotated-maf $annotated_maf"   if ( -s $annotated_maf );

# Run maf2maf.pl to annotate variants
system( $wrapper_cmd ) == 0  or die "\nERROR: Failed to run $vep_wrap_script!\nCommand: $wrapper_cmd\n";

# Replace previous annotated MAF
#`cp $output_maf $annotated_maf` if ( -s $output_maf );


__DATA__

=head1 NAME
 
 msk-vep-wrapper.pl - This script stores and re-uses VEP-annotated MAF files. It is designed specifically for MSKCC cbioportal
 
=head1 SYNOPSIS
 
 perl msk-vep-wrapper.pl --help
 perl msk-vep-wrapper.pl --input-maf test.maf --output-maf test.vep.maf --maf2maf /ssd-data/cmo/opt/vcf2maf/maf2maf.pl --tmp-dir /ssd-data/portal-tmp-dir
 
=head1 OPTIONS
 
 --input-maf      Path to input file in MAF format
 --output-maf     Path to output MAF file [Default: STDOUT]
 --tmp-dir        Folder to retain intermediate VCFs/MAFs after runtime [Default: usually under /tmp]
 --tum-depth-col  Name of MAF column for read depth in tumor BAM [t_depth]
 --tum-rad-col    Name of MAF column for reference allele depth in tumor BAM [t_ref_count]
 --tum-vad-col    Name of MAF column for variant allele depth in tumor BAM [t_alt_count]
 --nrm-depth-col  Name of MAF column for read depth in normal BAM [n_depth]
 --nrm-rad-col    Name of MAF column for reference allele depth in normal BAM [n_ref_count]
 --nrm-vad-col    Name of MAF column for variant allele depth in normal BAM [n_alt_count]
 --retain-cols    Comma-delimited list of columns to retain from the input MAF
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing variant_effect_predictor.pl
 --vep-data       VEP's base cache/plugin directory
 --vep-forks      Number of forked processes to use when running VEP [4]
 --ref-fasta      Reference FASTA file
 --maf2maf        Path to script maf2maf.pl [/ssd-data/cmo/opt/vcf2maf/maf2maf.pl]
 --help           Print a brief help message and quit
 --man            Print the detailed manual
 
=head1 DESCRIPTION
 
 This script stores annotated MAF files, with which to run fast-vep-wrapper.pl to save time. It should be saved in the same directory as fast-vep-wrapper.pl.
 
=head1 AUTHORS
 
 Qingguo Wang
 
 Nikolaus Schultz Lab
 Memorial Sloan Kettering Cancer Center
 New York, NY 10065
 
=cut
