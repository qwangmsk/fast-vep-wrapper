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
my ( $ref_fasta, $config_file, $maf2maf );


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
'maf2maf=s' => \$maf2maf,
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );


( defined $input_maf ) or die "ERROR: input-maf is not defined!\n";


# Check presence of configuration file
if( $config_file ) {
    ( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
    # Use configuration file to initialize variables
    ReadConfigFile( );
}
else{
    $config_file = "$FindBin::Bin/config.txt";
    if ( -e $config_file ){
        ReadConfigFile( );
    }
    # elsif ( -e 'config.txt' || -e $ENV{"HOME"}.'/config.txt' ) {
    #    $config_file = (-e $ENV{"HOME"}.'/config.txt') ? $ENV{"HOME"}.'/config.txt' : 'config.txt';
    #    ReadConfigFile( );
    # }
}


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
} else {
    warn "WARNING: $var_count variants in input MAF will be annotated\n";
}
%mutations = ();


######################################## Run maf2maf ####################################################

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
    $maf2maf_cmd .= " --retain-cols \"$retain_cols\""   if ( $retain_cols );
    $maf2maf_cmd .= " --custom-enst $custom_enst_file"  if ( $custom_enst_file );
    $maf2maf_cmd .= " --vep-path $vep_path"             if ( $vep_path );
    $maf2maf_cmd .= " --vep-data $vep_data"             if ( $vep_data );
    $maf2maf_cmd .= " --vep-forks $vep_forks"           if ( $vep_forks );
    $maf2maf_cmd .= " --ref-fasta $ref_fasta"           if ( $ref_fasta );
    
    system( $maf2maf_cmd ) == 0  or die "\nERROR: Failed to run maf2maf!\nCommand: $maf2maf_cmd\n";
    `rm $input_maf.uniq`;
}


###################################### Process Input MAF ##################################################

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
        if( $tum_depth_col and defined $input_maf_col_idx{ lc( $tum_depth_col ) } ){
            push( @kept_cols, lc( $tum_depth_col ) ) ;
        }elsif( defined $input_maf_col_idx{ t_depth } ){
            push( @kept_cols, 't_depth' ) ;
        }
        if( $nrm_depth_col and defined $input_maf_col_idx{lc( $nrm_depth_col )} ){
            push ( @kept_cols, lc( $nrm_depth_col ) ) ;
        }elsif( defined $input_maf_col_idx{ n_depth } ){
            push ( @kept_cols, 'n_depth' ) ;
        }
        if( $tum_rad_col   and defined $input_maf_col_idx{lc( $tum_rad_col )} ){
            push ( @kept_cols, lc( $tum_rad_col ) ) ;
        }elsif( defined $input_maf_col_idx{ t_ref_count } ){
            push ( @kept_cols, 't_ref_count' ) ;
        }
        if( $tum_vad_col and defined $input_maf_col_idx{lc( $tum_vad_col )} ){
            push ( @kept_cols, lc( $tum_vad_col ) );
        }elsif( defined $input_maf_col_idx{ t_alt_count } ){
            push ( @kept_cols, 't_alt_count' );
        }
        if( $nrm_rad_col and defined $input_maf_col_idx{lc( $nrm_rad_col )} ){
            push ( @kept_cols, lc( $nrm_rad_col ) );
        }elsif( defined $input_maf_col_idx{ n_ref_count } ){
            push ( @kept_cols, 'n_ref_count' );
        }
        if( $nrm_vad_col and defined $input_maf_col_idx{lc( $nrm_vad_col )} ){
            push ( @kept_cols, lc( $nrm_vad_col ) );
        }elsif( defined $input_maf_col_idx{ n_alt_count } ){
            push ( @kept_cols, 'n_alt_count' );
        }
        if( $retain_cols ){
            map{ my $c = lc; push( @kept_cols, $c ) if ( !exists $force_new_cols{ $c } ) } split( ",", $retain_cols );
        }
    }
    else {
        my ( $chr, $pos, $ref, $al1, $al2, $t_id ) = map{ my $c = lc; ( defined $input_maf_col_idx{ $c } ? $cols[ $input_maf_col_idx{ $c } ] : undef )} qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode );
        
        # Make sure we have at least 1 variant allele. If 1 is unset, set it to the reference allele
        next if( !$al1 and !$al2 );
        
        $al1 = $ref unless( $al1 );
        $al2 = $ref unless( $al2 );
        # To simplify setting tumor genotype later, ensure that $al2 is always non-REF
        ( $al1, $al2 ) = ( $al2, $al1 ) if( $al2 eq $ref );
        # Create a key for this variant using Chromosome:Start_Position:Tumor_Sample_Barcode:Reference_Allele:Variant_Allele
        my $key = join( ":", ( $chr, $pos, $t_id, $ref, $al1, $al2 ) );
        # Count duplicate mutations
        foreach( 0..100 ) {
           if( !exists $input_maf_data{ $key."\t$_" } ) {
               $key .= "\t$_";
               last;
           }
        }
        
        # Store values for this variant into a hash, adding column names to the key
        foreach my $c ( @kept_cols ) {
            $input_maf_data{ $key }{ $c } = "";
            $input_maf_data{ $key }{ $c } = $cols[ $input_maf_col_idx{ $c } ] if( defined $input_maf_col_idx{ $c } && defined $cols[ $input_maf_col_idx{ $c } ] );
        }
        $input_maf_data{ $key }{ matched_norm_sample_barcode } = "NORMAL" unless ( $input_maf_data{ $key }{ matched_norm_sample_barcode } );
        
        $key = join( ":", ( $chr, $pos, $ref, $al1, $al2 ) );
        if ( defined $input_maf_data{ $key }{ '(AllTumorBarcodeTogether)+' } ) {
            $input_maf_data{ $key }{ '(AllTumorBarcodeTogether)+' } .= ','.$cols[ $input_maf_col_idx{ tumor_sample_barcode } ];
        } else {
            $input_maf_data{ $key }{ '(AllTumorBarcodeTogether)+' }  = $cols[ $input_maf_col_idx{ tumor_sample_barcode } ];
        }
        
    }
}
$input_maf_fh->close;

# Fetch column header
my $maf_header;
if ( $annotated_maf ) {
    $maf_header = `grep ^Hugo_Symbol $annotated_maf`;
    chomp $maf_header;
}
    
# Sanity check
if ( ! $maf_header ) {
    $maf_header = `grep ^Hugo_Symbol $input_maf.uniq.maf`;
    chomp $maf_header ;
}
elsif ( -e "$input_maf.uniq.maf" ) {
    my $new_maf_header = `grep ^Hugo_Symbol $input_maf.uniq.maf`;
    chomp $new_maf_header ;
    die "ERROR: The format of new annotation is not consistent with $annotated_maf\nAnnotated MAF:\n$maf_header\nNew annotation:\n$new_maf_header\n" if ( $maf_header ne $new_maf_header );
}


my $output_maf_fh = *STDOUT;
if ( $output_maf ) {
    $output_maf_fh = IO::File->new( $output_maf, ">" ) or die "ERROR: Couldn't open file: $output_maf\n";
}

PropagateAnnotation ( "$annotated_maf", 1 ) if ( $annotated_maf );
if ( -e "$input_maf.uniq.maf" ){
    PropagateAnnotation ( "$input_maf.uniq.maf", ( $annotated_maf ? 0 : 1 ) );
    `rm $input_maf.uniq.maf` if ( -e "$input_maf.uniq.maf" );
}
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
    
        my ( $chr, $pos, $ref, $al1, $al2 ) = map{ my $c = lc; ( defined $output_maf_col_idx{ $c } ? $cols[ $output_maf_col_idx{ $c } ] : undef )} qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 );
    
        # Make sure we have at least 1 variant allele. If 1 is unset, set it to the reference allele
        next if( !$al1 and !$al2 );
    
        $al1 = $ref unless( $al1 );
        $al2 = $ref unless( $al2 );
        # To simplify setting tumor genotype later, ensure that $al2 is always non-REF
        ( $al1, $al2 ) = ( $al2, $al1 ) if( $al2 eq $ref );
        # Create a key for this variant using Chromosome:Start_Position:Tumor_Sample_Barcode:Reference_Allele:Variant_Allele
        my $key = join( ":", ( $chr, $pos, $ref, $al1, $al2 ) );
        next if ( !defined $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} );
    
        my @t_ids = split( ",", $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} );
        $input_maf_data{$key}{'(AllTumorBarcodeTogether)+'} = "";
        my %flag = ();
        foreach my $t ( @t_ids ) {
            next if ( exists $flag{ $t } );
            $flag{ $t } = 1;
            $cols[ $output_maf_col_idx{ tumor_sample_barcode } ] = $t;
            foreach( 0..100 ){
                $key = join( ":", ( $chr, $pos, $t, $ref, $al1, $al2 ) ) . "\t$_";
                last if ( !exists $input_maf_data{ $key });
                
                foreach my $c ( @kept_cols ){
                    $cols[ $output_maf_col_idx{ $c } ] = $input_maf_data{ $key }{ $c };
                }
                
                if( defined $tum_rad_col and exists $input_maf_data{ $key }{ lc( $tum_rad_col ) } ){
                    $cols[ $output_maf_col_idx{ t_ref_count } ] = $input_maf_data{ $key }{ lc( $tum_rad_col ) } ;
                }elsif( exists $input_maf_data{ $key }{ t_ref_count }  ){
                    $cols[ $output_maf_col_idx{ t_ref_count } ] = $input_maf_data{ $key }{ t_ref_count } ;
                }
                if( defined $tum_vad_col and exists $input_maf_data{ $key }{ lc( $tum_vad_col ) } ){
                    $cols[ $output_maf_col_idx{ t_alt_count } ] = $input_maf_data{ $key }{ lc( $tum_vad_col ) } ;
                }elsif( exists $input_maf_data{ $key }{ t_alt_count } ){
                    $cols[ $output_maf_col_idx{ t_alt_count } ] = $input_maf_data{ $key }{ t_alt_count } ;
                }
                if( defined $nrm_rad_col and exists $input_maf_data{ $key }{ lc( $nrm_rad_col ) } ){
                    $cols[ $output_maf_col_idx{ n_ref_count } ] = $input_maf_data{ $key }{ lc( $nrm_rad_col ) } ;
                }elsif( exists $input_maf_data{ $key }{ n_ref_count } ){
                    $cols[ $output_maf_col_idx{ n_ref_count } ] = $input_maf_data{ $key }{ n_ref_count } ;
                }
                if( defined $nrm_vad_col and exists $input_maf_data{ $key }{ lc( $nrm_vad_col ) } ){
                    $cols[ $output_maf_col_idx{ n_alt_count } ] = $input_maf_data{ $key }{ lc( $nrm_vad_col ) } ;
                }elsif( exists $input_maf_data{ $key }{ n_alt_count } ){
                    $cols[ $output_maf_col_idx{ n_alt_count } ] = $input_maf_data{ $key }{ n_alt_count } ;
                }
                
                if( !defined $tum_depth_col or !exists $input_maf_col_idx{lc( $tum_depth_col )} or !defined $cols[ $output_maf_col_idx{ t_depth } ] or
                (( $tum_depth_col and exists $input_maf_col_idx{lc( $tum_depth_col )} and defined $cols[ $output_maf_col_idx{ t_depth } ] and
                (( $cols[ $output_maf_col_idx{ t_ref_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ t_ref_count }] > $cols[ $output_maf_col_idx{ t_depth } ] ) or
                ( $cols[ $output_maf_col_idx{ t_alt_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ t_alt_count }] > $cols[ $output_maf_col_idx{ t_depth } ] ) or
                ( $cols[ $output_maf_col_idx{ t_ref_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ t_alt_count }] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ t_ref_count }] + $cols[$output_maf_col_idx{ t_alt_count }] > $cols[ $output_maf_col_idx{ t_depth } ] ))))) {
                    $cols[ $output_maf_col_idx{ t_depth } ]  = 0;
                    $cols[ $output_maf_col_idx{ t_depth } ] += $cols[$output_maf_col_idx{ t_ref_count }] if ( $cols[$output_maf_col_idx{ t_ref_count }] =~ m/^\d+$/ );
                    $cols[ $output_maf_col_idx{ t_depth } ] += $cols[$output_maf_col_idx{ t_alt_count }] if ( $cols[$output_maf_col_idx{ t_alt_count }] =~ m/^\d+$/ );
                }else{
                    if( defined $tum_depth_col and defined $input_maf_data{$key}{ lc( $tum_depth_col ) } ){
                        $cols[ $output_maf_col_idx{ t_depth } ] =  $input_maf_data{ $key }{ lc( $tum_depth_col ) };
                    }elsif ( exists $input_maf_data{ $key }{ t_depth } ){
                        $cols[ $output_maf_col_idx{ t_depth } ] =  $input_maf_data{ $key }{ t_depth };
                    }else{
                        $cols[ $output_maf_col_idx{ t_depth } ] =  "";
                    }
                }
                
                if( !defined $nrm_depth_col or !exists $input_maf_col_idx{lc( $nrm_depth_col )} or !defined $cols[ $output_maf_col_idx{ n_depth } ] or
                (( $nrm_depth_col and exists $input_maf_col_idx{lc( $nrm_depth_col )} and defined $cols[ $output_maf_col_idx{ n_depth } ] and
                (( $cols[ $output_maf_col_idx{ n_ref_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ n_ref_count }] > $cols[ $output_maf_col_idx{ n_depth } ] ) or
                ( $cols[ $output_maf_col_idx{ n_alt_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ n_alt_count }] > $cols[ $output_maf_col_idx{ n_depth } ] ) or
                ( $cols[ $output_maf_col_idx{ n_ref_count } ] =~ m/^\d+$/ and $cols[$output_maf_col_idx{ n_alt_count }] =~ m/^\d+$/ and $cols[ $output_maf_col_idx{ n_ref_count } ] + $cols[$output_maf_col_idx{ n_alt_count }] > $cols[ $output_maf_col_idx{ n_depth } ] ))))) {
                    $cols[ $output_maf_col_idx{ n_depth } ]  = 0;
                    $cols[ $output_maf_col_idx{ n_depth } ] += $cols[$output_maf_col_idx{ n_ref_count }] if ( $cols[$output_maf_col_idx{ n_ref_count }] =~ m/^\d+$/ );
                    $cols[ $output_maf_col_idx{ n_depth } ] += $cols[$output_maf_col_idx{ n_alt_count }] if ( $cols[$output_maf_col_idx{ n_alt_count }] =~ m/^\d+$/ );
                }else{
                    if( defined $nrm_depth_col and defined $input_maf_data{$key} {lc( $nrm_depth_col ) } ){
                        $cols[ $output_maf_col_idx{ n_depth } ] =  $input_maf_data{ $key }{ lc( $nrm_depth_col ) };
                    }elsif( $input_maf_data{ $key }{ n_depth } ){
                        $cols[ $output_maf_col_idx{ n_depth } ] = $input_maf_data{ $key }{ n_depth };
                    }else{
                        $cols[ $output_maf_col_idx{ n_depth } ] =  "";
                    }
                }
                $output_maf_fh->print( join( "\t", @cols ) . "\n" );
            }
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
        my ( $chr, $pos, $ref, $al1, $al2, $t_id ) = map{ my $c = lc; ( defined $col_idx{$c} ? $cols[ $col_idx{$c} ] : undef )} qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode);
        
        # Make sure we have at least 1 variant allele. If 1 is unset, set it to the reference allele
        if( !$al1 and !$al2 ) {
            warn "WARNING: Skipping variant at $chr:$pos without any variant alleles specified!\n";
            next;
        }
        
        # If normal alleles are unset in the MAF (quite common), assume homozygous reference
        $al1 = $ref unless( $al1 );
        $al2 = $ref unless( $al2 );
        # To simplify setting tumor genotype later, ensure that $al2 is always non-REF
        ( $al1, $al2 ) = ( $al2, $al1 ) if( $al2 eq $ref );
        
        # Construct a hash key to filter out duplicate mutations.
        my $key = join( ":", ( $chr, $pos, $ref, $al1, $al2 ) );
        
        if ( ! exists $mutations { $key } ) {
            $mutations { $key } = 1;
            
            $cols[$t_idx] = 'TUMOR';
            $cols[$n_idx] = 'NORMAL' if ( defined $n_idx );
            
            $out_maf_fh->print( join( "\t", @cols) . "\n" ) if ( defined $out_maf_fh );
            $line_count++;
        }
    }
    $in_maf_fh->close;
    
    return $line_count;
}


################################### Read configuration file #########################################


sub ReadConfigFile {
    my %config;
    map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;
    
    # Use the configuration file to initialize variables
    $custom_enst_file = $config{ custom_enst_file } if ( !$custom_enst_file && exists $config{ custom_enst_file } );
    $vep_path         = $config{ vep_path }         if ( !$vep_path         && exists $config{ vep_path } );
    $vep_data         = $config{ vep_data }         if ( !$vep_data         && exists $config{ vep_data } );
    $ref_fasta        = $config{ ref_fasta }        if ( !$ref_fasta        && exists $config{ ref_fasta } );
    $vep_forks        = $config{ vep_forks }        if ( !$vep_forks        && exists $config{ vep_forks } );
    $tmp_dir          = $config{ tmp_dir }          if ( !$tmp_dir          && exists $config{ tmp_dir } );
    $maf2maf          = $config{ vcf2maf_script }   if ( !$maf2maf          && exists $config{ vcf2maf_script } );
}


__DATA__

=head1 NAME
 
 fast-vep-wrapper.pl - Quick annotation of variants in MAF files.
 
=head1 SYNOPSIS
 
 perl fast-vep-wrapper.pl --help
 perl fast-vep-wrapper.pl --input-maf test.maf --output-maf test.vep.maf --maf2maf /ssd-data/cmo/opt/vcf2maf/maf2maf.pl
 perl fast-vep-wrapper.pl --input-maf test.maf --output-maf test.vep.maf --config-file /home/someone/bin/config.txt
 
=head1 OPTIONS
 
 --input-maf      Path to input file in MAF format
 --output-maf     Path to output MAF file [Default: STDOUT]
 --annotated-maf  MAF file annotated previously using VEP
 --tmp-dir        Folder to retain intermediate VCFs/MAFs after runtime [Default: usually under /tmp]
 --tum-depth-col  Name of MAF column for read depth in tumor BAM [t_depth]
 --tum-rad-col    Name of MAF column for reference allele depth in tumor BAM [t_ref_count]
 --tum-vad-col    Name of MAF column for variant allele depth in tumor BAM [t_alt_count]
 --nrm-depth-col  Name of MAF column for read depth in normal BAM [n_depth]
 --nrm-rad-col    Name of MAF column for reference allele depth in normal BAM [n_ref_count]
 --nrm-vad-col    Name of MAF column for variant allele depth in normal BAM [n_alt_count]
 --retain-cols    Comma-delimited list of columns to retain from the input MAF
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing variant_effect_predictor.pl [/ssd-data/cmo/opt/vep/v79]
 --vep-data       VEP's base cache/plugin directory [/ssd-data/cmo/opt/vep/v79]
 --vep-forks      Number of forked processes to use when running VEP [4]
 --ref-fasta      Reference FASTA file [/ssd-data/cmo/opt/vep/v79/homo_sapiens/79_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa]
 --config-file    A configuration file to store paths of vep, ref_fasta, and maf2maf, in addition to other arguments
 --maf2maf        Path to script maf2maf.pl [/ssd-data/cmo/opt/vcf2maf/maf2maf.pl]
 --help           Print a brief help message and quit
 --man            Print the detailed manual
 
=head1 DESCRIPTION
 
 This script speeds up MAF file annotation using VEP.
 
=head1 AUTHORS
 
 Qingguo Wang
 
 Nikolaus Schultz Lab
 Memorial Sloan Kettering Cancer Center
 New York, NY 10065

=head1 ACKNOWLEDGEMENTS
 
 Thank Sumit Middha, Cyriac Kandoth, and Frederick Criscuolo for helpful suggestions and discussion
 
=cut
