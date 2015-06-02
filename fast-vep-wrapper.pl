#!/usr/bin/env perl

# fast-vep-wrapper - Quick annotation of variants in a MAF file


use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use FindBin;
use lib "$FindBin::Bin";
use Cwd;


# Columns that should never be overridden since they are results of re-annotation
my %force_new_cols = map{ my $c = lc; ( $c, 1 )} qw( Hugo_Symbol Entrez_Gene_Id NCBI_Build
Chromosome Start_Position End_Position Strand Variant_Classification Variant_Type
Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode
Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1
Match_Norm_Validation_Allele2 HGVSc HGVSp HGVSp_Short Transcript_ID Exon_Number t_depth
t_ref_count t_alt_count n_depth n_ref_count n_alt_count all_effects Allele Gene Feature
Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons
Existing_variation ALLELE_NUM DISTANCE STRAND SYMBOL SYMBOL_SOURCE HGNC_ID BIOTYPE CANONICAL
CCDS ENSP SWISSPROT TREMBL UNIPARC RefSeq SIFT PolyPhen EXON INTRON DOMAINS GMAF AFR_MAF
AMR_MAF ASN_MAF EUR_MAF AA_MAF EA_MAF CLIN_SIG SOMATIC PUBMED MOTIF_NAME MOTIF_POS
HIGH_INF_POS MOTIF_SCORE_CHANGE );

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

my ( $man, $help ) = ( 0, 0 );
my ( $input_maf, $output_maf, $annotated_maf, $tmp_dir, $custom_enst_file );
my ( $tum_depth_col, $tum_rad_col, $tum_vad_col );
my ( $nrm_depth_col, $nrm_rad_col, $nrm_vad_col );
my ( $retain_cols, $vep_path, $vep_data, $vep_forks );
my ( $ref_fasta, $config_file );


GetOptions(
'help!' => \$help,
'man!' => \$man,
'input-maf=s' => \$input_maf,
'output-maf=s' => \$output_maf,
'annotated-maf=s' => \$annotated_maf,
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
'config-file=s' => \$config_file,
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );


( defined $input_maf ) or die "ERROR: input-maf is defined!\n";


################################# Filter out duplicate variants ##########################################

my %mutations;

GetUniqVariants( $annotated_maf )  if ( $annotated_maf );

# Extract unique variants from input maf file
my $maf_fh = IO::File->new( "$input_maf.uniq", ">" ) or die "ERROR: Couldn't create temporary file $input_maf.uniq\n";
my $var_count = GetUniqVariants( $input_maf, $maf_fh );
$maf_fh->close;

if( $var_count == 0 ) {
    warn "WARNING: 0 unannotated variants in input MAF\n";
    `rm "$input_maf.uniq"`;
}
%mutations = ();

######################################## Run maf2maf ####################################################

# Check configuration file
if ($config_file) {
    ( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
}else{
    $config_file = "$FindBin::Bin/config.txt";
    if (!-e $config_file){
        die "ERROR: Could not find configuration file config.txt in home/working/program directory\n" if (!-e 'config.txt' && !-e $ENV{"HOME"}.'/config.txt');
        $config_file = (-e $ENV{"HOME"}.'/config.txt') ? $ENV{"HOME"}.'/config.txt' : 'config.txt'
    }
}

# Read configuration file
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;

my $maf2maf= $config{ 'vcf2maf_script' };
$custom_enst_file   = $config{ 'custom_enst_file' }     if ( !$custom_enst_file  && defined $config{ 'custom_enst_file' } );
$vep_path           = $config{ 'vep_path' }             if ( !$vep_path  && defined $config{ 'vep_path' } );
$vep_data           = $config{ 'vep_data' }             if ( !$vep_data  && defined $config{ 'vep_data' } );
$ref_fasta          = $config{ 'ref_fasta' }            if ( !$ref_fasta && defined $config{ 'ref_fasta' } );
$vep_forks          = $config{ 'vep_forks' }            if ( !$vep_forks && defined $config{ 'vep_forks' } );
$tmp_dir            = $config{ 'tmp_dir' }              if ( !$tmp_dir   && defined $config{ 'tmp_dir' } );

# Run maf2maf.pl to annotate variants
if ( $var_count > 0 ) {
    
    my $maf2maf_cmd = "perl $maf2maf --input-maf $input_maf.uniq --output-maf $input_maf.uniq.maf";
    
    $maf2maf_cmd .= " --tmp-dir $tmp_dir"               if ( $tmp_dir );
    $maf2maf_cmd .= " --tum-depth-col $tum_depth_col"   if ( $tum_depth_col );
    $maf2maf_cmd .= " --tum-rad-col $tum_rad_col"       if ( $tum_rad_col );
    $maf2maf_cmd .= " --tum-vad-col $tum_vad_col"       if ( $tum_vad_col );
    $maf2maf_cmd .= " --nrm-depth-col $nrm_depth_col"   if ( $nrm_depth_col );
    $maf2maf_cmd .= " --nrm-rad-col $nrm_rad_col"       if ( $nrm_rad_col );
    $maf2maf_cmd .= " --nrm-vad-col $nrm_vad_col"       if ( $nrm_vad_col );
    $maf2maf_cmd .= " --retain-cols $retain_cols"       if ( $retain_cols );
    $maf2maf_cmd .= " --custom-enst $custom_enst_file"  if ( $custom_enst_file );
    $maf2maf_cmd .= " --vep-path $vep_path"             if ( $vep_path );
    $maf2maf_cmd .= " --vep-data $vep_data"             if ( $vep_data );
    $maf2maf_cmd .= " --vep-forks $vep_forks"           if ( $vep_forks );
    $maf2maf_cmd .= " --ref-fasta $ref_fasta"           if ( $ref_fasta );
    
    system( $maf2maf_cmd ) == 0  or die "\nERROR: Failed to run maf2maf!\nCommand: $maf2maf_cmd\n";
}

###################################### Assign Annotation ##################################################

# Parse the input MAF and fetch the data for columns that we need to retain/override
my $input_maf_fh = IO::File->new( $input_maf ) or die "ERROR: Couldn't open file: $input_maf\n";

my %input_maf_col_idx = (); # Hash to map column names to column indexes
my %input_maf_data = (); # fetch and override columns from the input MAF that user wants to retain
my @kept_cols;

while( my $line = $input_maf_fh->getline ) {
    
    next if( $line =~ m/^#/ ); # Skip comments
        
    # Do a thorough removal of carriage returns, line feeds, prefixed/suffixed whitespace
    my @cols = map{s/^\s+|\s+$|\r|\n//g; $_} split( /\t/, $line );
    
    # Parse the header line to map column names to their indexes
    if( $line =~ m/^(Hugo_Symbol|Chromosome)/ ) {
        my $idx = 0;
        map{ my $c = lc; $input_maf_col_idx{$c} = $idx; ++$idx } @cols;
        
        @kept_cols= qw( tumor_sample_barcode matched_norm_sample_barcode );
        push (@kept_cols, lc( $tum_depth_col )) if( $tum_depth_col && defined $input_maf_col_idx{lc $tum_depth_col} );
        push (@kept_cols, lc( $nrm_depth_col )) if( $nrm_depth_col && defined $input_maf_col_idx{lc $nrm_depth_col} );
        push (@kept_cols, lc( $tum_rad_col ))   if( $tum_rad_col   && defined $input_maf_col_idx{lc $tum_rad_col} );
        push (@kept_cols, lc( $tum_vad_col ))   if( $tum_vad_col   && defined $input_maf_col_idx{lc $tum_vad_col} );
        push (@kept_cols, lc( $nrm_rad_col ))   if( $nrm_rad_col   && defined $input_maf_col_idx{lc $nrm_rad_col} );
        push (@kept_cols, lc( $nrm_vad_col ))   if( $nrm_vad_col   && defined $input_maf_col_idx{lc $nrm_vad_col} );
         
        if( $retain_cols ){
            map{my $c = lc; push( @kept_cols, $c ) if ( !$force_new_cols{ $c } ) } split( ",", $retain_cols );
        }
    }
    else {
        # Figure out which of the tumor alleles is non-reference
        my ( $ref, $al1, $al2 ) = map{ my $c = lc; ( defined $input_maf_col_idx{$c} ? $cols[$input_maf_col_idx{$c}] : "" ) } qw( Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 );
        my $var_allele = (( defined $al1 and $al1 and $al1 ne $ref ) ? $al1 : $al2 );
        
        # Create a key for this variant using Chromosome:Start_Position:Tumor_Sample_Barcode:Reference_Allele:Variant_Allele
        my $key = join( ":", ( map{ my $c = lc; $cols[$input_maf_col_idx{$c}] } qw( Chromosome Start_Position Tumor_Sample_Barcode Reference_Allele )), $var_allele );
        
        # Store values for this variant into a hash, adding column names to the key
        foreach my $c ( @kept_cols ) {
            $input_maf_data{$key}{$c} = "";
            if( defined $input_maf_col_idx{$c} and defined $cols[$input_maf_col_idx{$c}] ) {
                $input_maf_data{$key}{$c} = $cols[$input_maf_col_idx{$c}];
            }
        }
        
        $key = join( ":", ( map{ my $c = lc; $cols[$input_maf_col_idx{$c}] } qw( Chromosome Start_Position Reference_Allele )), $var_allele );
        if ( defined $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} ) {
            $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} .= ','.$cols[ $input_maf_col_idx{ tumor_sample_barcode } ];
        } else {
            $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'}  = $cols[ $input_maf_col_idx{ tumor_sample_barcode } ];
        }
        
    }
}
$input_maf_fh->close;

# Fetch column header
my $maf_header;
$maf_header = `grep ^Hugo_Symbol $annotated_maf` if ( $annotated_maf );
chomp $maf_header ;

# Sanity check
if ( ! $maf_header ) {
    $maf_header = `grep ^Hugo_Symbol $input_maf.uniq.maf`;
    chomp $maf_header ;
}
elsif ( -e "$input_maf.uniq.maf" ) {
    my $new_maf_header = `grep ^Hugo_Symbol $input_maf.uniq.maf`;
    chomp $new_maf_header ;
    die "ERROR: The format of new VEP annotation is inconsistent with $annotated_maf\n" if ( $maf_header ne $new_maf_header );
}


my $output_maf_fh = *STDOUT;
if ( $output_maf ) {
    $output_maf_fh = IO::File->new( $output_maf, ">" ) or die "ERROR: Couldn't open file: $output_maf\n";
}

PropagateAnnotation ( "$annotated_maf", 1 ) if ( $annotated_maf );
PropagateAnnotation ( "$input_maf.uniq.maf", ( $annotated_maf ? 0 : 1 ) ) if ( -e "$input_maf.uniq.maf" );

$output_maf_fh->close;

########################## Read VEP annotation to be assign to variants in input MAF ##############################


sub PropagateAnnotation {
    
    my $anno_maf  = shift;
    my $print_header = shift;
    
    # Retain/override data in MAF
    my $anno_maf_fh = IO::File->new( $anno_maf ) or die "ERROR: Couldn't open file: $anno_maf\n";
    my %output_maf_col_idx = (); # Hash to map column names to column indexes
    
    while( my $line = $anno_maf_fh->getline ) {
        
        # Do a thorough removal of carriage returns, line feeds, prefixed/suffixed whitespace
        my @cols = map{ s/^\s+|\s+$|\r|\n//g; $_ } split( /\t/, $line );
        
        # Copy comment lines to the new MAF unchanged
        if( $line =~ m/^#/ ) {
            $output_maf_fh->print( $line ) if ( $print_header );
            next;
        }
        # Print the MAF header prepared earlier, but also create a hash with column indexes
        elsif( $line =~ m/^Hugo_Symbol/ ) {
            my $idx = 0;
            map{ my $c = lc; $output_maf_col_idx{$c} = $idx; ++$idx } ( @cols );
            $output_maf_fh->print( "$maf_header\n" ) if ( $print_header );
            next;
        }
        
        # For all other lines, insert the data collected from the original input MAF
        my $key = join( ":", map{ my $c = lc; $cols[$output_maf_col_idx{$c}] } qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele2 ));
        next if ( $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} eq "");
    
        my @t_ids = split( ",", $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} );
        $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} = "";

        foreach my $t ( @t_ids ) {
            $cols[$output_maf_col_idx{tumor_sample_barcode}] = $t;
            $key = join( ":", map{ my $c = lc; $cols[$output_maf_col_idx{$c}] } qw( Chromosome Start_Position Tumor_Sample_Barcode Reference_Allele Tumor_Seq_Allele2));
            foreach my $c ( @kept_cols ){
                $cols[$output_maf_col_idx{$c}] = $input_maf_data{$key}{$c} if( defined $input_maf_data{$key}{$c} );
            }
            
            $cols[ $output_maf_col_idx{lc('t_ref_count')} ] = ( defined $tum_rad_col && defined $input_maf_data{$key}{$tum_rad_col} ) ? $input_maf_data{$key}{$tum_rad_col} : "";
            $cols[ $output_maf_col_idx{lc('t_alt_count')} ] = ( defined $tum_vad_col && defined $input_maf_data{$key}{$tum_vad_col} ) ? $input_maf_data{$key}{$tum_vad_col} : "";
            $cols[ $output_maf_col_idx{lc('n_ref_count')} ] = ( defined $nrm_rad_col && defined $input_maf_data{$key}{$nrm_rad_col} ) ? $input_maf_data{$key}{$nrm_rad_col} : "";
            $cols[ $output_maf_col_idx{lc('n_alt_count')} ] = ( defined $nrm_vad_col && defined $input_maf_data{$key}{$nrm_vad_col} ) ? $input_maf_data{$key}{$nrm_vad_col} : "";
            
            $cols[ $output_maf_col_idx{lc('t_depth')} ] = ( $cols[$output_maf_col_idx{lc('t_ref_count')}] =~ /^\d+$/ && $cols[$output_maf_col_idx{lc('t_alt_count')}] =~ /^\d+$/ ) ? ( $cols[$output_maf_col_idx{lc('t_ref_count')}] + $cols[$output_maf_col_idx{lc('t_alt_count')}] ) : "";
            $cols[ $output_maf_col_idx{lc('n_depth')} ] = ( $cols[$output_maf_col_idx{lc('n_ref_count')}] =~ /^\d+$/ && $cols[$output_maf_col_idx{lc('n_alt_count')}] =~ /^\d+$/ ) ? ( $cols[$output_maf_col_idx{lc('n_ref_count')}] + $cols[$output_maf_col_idx{lc('n_alt_count')}] ) : "";
            
            $output_maf_fh->print( join( "\t", @cols ) . "\n" );
        }
    }
    
    $anno_maf_fh->close;
    
}


###################################### Get unique variants ##################################################


sub GetUniqVariants {
    
    my $in_maf     = shift;
    my $out_maf_fh = shift;
    
    my $line_count = 0;
    my %col_idx = (); # Hash to map column names to column indexes
    my ($t_idx, $n_idx);
    
    # Parse through each variant in the MAF
    my $in_maf_fh = IO::File->new( $in_maf ) or die "ERROR: Couldn't open file: $in_maf\n";
    
    while( my $line = $in_maf_fh->getline ) {
        
        # Skip comment lines
        next if( $line =~ m/^#/ );
            
        # Instead of a chomp, do a thorough removal of carriage returns, line feeds, and prefixed/suffixed whitespace
        my @cols = map{s/^\s+|\s+$|\r|\n//g; $_} split( /\t/, $line );
        
        # Parse the header line to map column names to their indexes
        if( $line =~ m/^(Hugo_Symbol|Chromosome)/ ) {
            my $idx = 0;
            
            # Fetch the column names and do some sanity checks (don't be case-sensitive)
            map{ my $c = lc; $col_idx{$c} = $idx; ++$idx; } @cols;
            map{ my $c = lc; ( defined $col_idx{$c} ) or die "ERROR: $_ is a required MAF column!\n" } qw( Chromosome Start_Position Reference_Allele Tumor_Sample_Barcode );
            ( defined $col_idx{tumor_seq_allele1} or defined $col_idx{tumor_seq_allele2} ) or die "ERROR: At least one MAF column for Tumor_Seq_Allele must be defined!\n";
            
            # Fetch all tumor-normal paired IDs from the MAF, doing some whitespace cleanup in the same step
            $t_idx = $col_idx{tumor_sample_barcode};
            $n_idx = $col_idx{matched_norm_sample_barcode} if( defined $col_idx{matched_norm_sample_barcode} );
            
            $out_maf_fh->print( $line ) if ( defined $out_maf_fh );
            next;
        }
        
        # Print an error if we got to this point without parsing a header line, and increment a counter for all non-header lines
        ( %col_idx ) or die "ERROR: Couldn't find a header line in the MAF: $in_maf";
        
        # For a variant in the MAF, parse out the bare minimum data needed by a VCF
        my ( $chr, $pos, $ref, $al1, $al2, $t_id, $n_id, $n_al1, $n_al2 ) = map{ my $c = lc; ( defined $col_idx{$c} ? $cols[$col_idx{$c}] : undef )} qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 );
        
        # If normal alleles are unset in the MAF (quite common), assume homozygous reference
        $n_al1 = $ref unless( $n_al1 );
        $n_al2 = $ref unless( $n_al2 );
        
        # Make sure we have at least 1 variant allele. If 1 is unset, set it to the reference allele
        if( !$al1 and !$al2 ) {
            warn "WARNING: Skipping variant at $chr:$pos without any variant alleles specified!\n";
            next;
        }
        
        $al1 = $ref unless( $al1 );
        $al2 = $ref unless( $al2 );
        # To simplify setting tumor genotype later, ensure that $al2 is always non-REF
        ( $al1, $al2 ) = ( $al2, $al1 ) if( $al2 eq $ref );
        # Do the same for the normal alleles, though it makes no difference if both are REF
        ( $n_al1, $n_al2 ) = ( $n_al2, $n_al1 ) if( $n_al2 eq $ref );
        
        # Construct a hash key to filter out duplicate mutations.
        my $hash_key = "$chr\t$pos\t$ref\t$al1\t$al2\t$n_al1\t$n_al2";
        
        if ( ! exists $mutations { $hash_key } ) {
            $mutations { $hash_key } = 1;
            
            $cols[$t_idx] = 'TUMOR';
            $cols[$n_idx] = 'NORMAL' if ( defined $n_idx );
            
            $out_maf_fh->print( join( "\t", @cols) . "\n" ) if ( defined $out_maf_fh );
            $line_count++;
        }
    }
    $in_maf_fh->close;
    
    return $line_count;
}


__DATA__

=head1 NAME
 
 fast-vep-wrapper.pl - Annotate the effects of variants in a MAF efficiently.
 
=head1 SYNOPSIS
 
 perl fast-vep-wrapper.pl --help
 perl fast-vep-wrapper.pl --input-maf test.maf --output-maf test.vep.maf  # Suppose configuration file config.txt is in predefined path
 perl fast-vep-wrapper.pl --input-maf test.maf --output-maf test.vep.maf --config-file /home/someone/bin/fast-vep-runner/config.txt
 
=head1 OPTIONS
 
 --input-maf      Path to input file in MAF format
 --output-maf     Path to output MAF file [Default: STDOUT]
 --annotated-maf  Existing MAF file annotated earlier using VEP
 --tmp-dir        Folder to retain intermediate VCFs/MAFs after runtime [Default: usually under /tmp]
 --tum-depth-col  Name of MAF column for read depth in tumor BAM [t_depth]
 --tum-rad-col    Name of MAF column for reference allele depth in tumor BAM [t_ref_count]
 --tum-vad-col    Name of MAF column for variant allele depth in tumor BAM [t_alt_count]
 --nrm-depth-col  Name of MAF column for read depth in normal BAM [n_depth]
 --nrm-rad-col    Name of MAF column for reference allele depth in normal BAM [n_ref_count]
 --nrm-vad-col    Name of MAF column for variant allele depth in normal BAM [n_alt_count]
 --retain-cols    Comma-delimited list of columns to retain from the input MAF [Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID]
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing variant_effect_predictor.pl [/ssd-data/cmo/opt/vep/v79]
 --vep-data       VEP's base cache/plugin directory [/ssd-data/cmo/opt/vep/v79]
 --vep-forks      Number of forked processes to use when running VEP [4]
 --ref-fasta      Reference FASTA file [/ssd-data/cmo/opt/vep/v79/homo_sapiens/79_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa]
 --config-file    A configuration file to store paths of vep, ref_fasta, and maf2maf. Default is config.txt in the same directory as this script. If not provided, will look for it in current directory and home directory.
 --help           Print a brief help message and quit
 --man            Print the detailed manual
 
=head1 DESCRIPTION
 
 This script speeds up MAF file annotation using VEP. It extracts unique variants from input MAF file. Here, 'unique' means unique combination of Chromosome, Start_Position, Reference_Allele and Tumor_Allele. If an annotated MAF file is provided, variants in the annotated MAF will not be re-annotated. Only new variants in the MAF file are processed (using vcf2maf: https://github.com/ckandoth/vcf2maf). For new variants, 'TUMOR' and 'NORMAL' are used as'Tumor_Sample_Barcode' and 'Matched_Norm_Sample_Barcode' to reduce multithreading overhead.
 Derived from vcf2maf (https://github.com/ckandoth/vcf2maf), this script aims to improve annotation efficiency that is lacking in vcf2maf. By reusing previous annotation and improving parallel efficiency, this script reduces time for annotating large MAF files from hours to minutes.
 
=head1 AUTHORS
 
 Qingguo Wang (josephw10000@gmail.com)

=head1 ACKNOWLEDGEMENTS
 
 Thank Sumit Middha, Frederick Criscuolo, and Cyriac Kandoth for suggestion and discussion
 
=cut