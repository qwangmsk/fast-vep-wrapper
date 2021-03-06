#!/usr/bin/perl -w

# run-vep-wrapper - An application-level wrapper script for annotating MAF files using VEP.


use strict;
use warnings FATAL => 'all';
use File::Basename;
use Getopt::Std;
use Getopt::Long;
use Pod::Usage qw( pod2usage );
use FindBin;
use lib "$FindBin::Bin";
use Cwd;


################################ Read input parameters ###############################

my ( $vep_wrap_script );
my ( $vep_path, $vep_data, $vep_forks );
my ( $ref_fasta, $custom_enst_file );
my ( $input, $input_filename, $output_maf, $annotated_maf );
my ( $tmp_dir, $depth_col_file);
my ( $config_file, $use_cluster );
my ( $help, $man );


GetOptions
(
'i|input=s'         => \$input,
'o|output-maf=s'    => \$output_maf,
'a|annotated-maf=s' => \$annotated_maf,
't|tmp-dir=s'       => \$tmp_dir,
'u|use-cluster=s'   => \$use_cluster,
'v|vep-forks=s'     => \$vep_forks,
'c|config-file=s'   => \$config_file,
'h|help|?'          => \$help,
'man!'              => \$man,
)or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );


# Check input and output arguments
( defined $input ) or die "ERROR: Missing input argument\n";
( -e $input ) or die "Error: $input does not exist\n";
if ( $output_maf ){
    die "ERROR: Output is a directory\n" if (-d $output_maf);
}

# Check configuration file
if ( defined $config_file ) {
    ( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
}else{
    $config_file = "$FindBin::Bin/config.txt";
    ( -e $config_file ) or die "ERROR: Could not find configuration file config.txt\n";
}


####################### Read configuration file ########################

my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;

$vep_path         = $config{ vep_path };
$vep_data         = $config{ vep_data };
$ref_fasta        = $config{ ref_fasta };
$custom_enst_file = $config{ custom_enst_file }  if ( defined $config{ custom_enst_file } );
$vep_forks        = $config{ vep_forks }         if ( !defined $vep_forks );
$tmp_dir          = $config{ tmp_dir }           if ( !defined $tmp_dir && defined $config{ tmp_dir } );
$vep_wrap_script  = $config{ vep_wrap_script };
$depth_col_file   = $config{ depth_def_file };
$input_filename   = $config{ input_filename };
$output_maf       = $config{ output_filename }   if ( !defined $output_maf    && defined $config{ output_filename } );
$annotated_maf    = $config{ annotated_maf }     if ( !defined $annotated_maf && defined $config{ annotated_maf } );
$use_cluster      = $config{ use_cluster }       if ( !defined $use_cluster   && defined $config{ use_cluster } );


# Check VEP wrapper scripts
( defined $vep_wrap_script and -e $vep_wrap_script ) or die "Error: Configuration file does not define vep_wrap_script\n";

# Check TMP dir
`mkdir -p $tmp_dir` if( defined $tmp_dir && !-e $tmp_dir );

# Check file existance
( -s $depth_col_file ) or die "Error: $depth_col_file does not exist\n";
# Get depth-related column names
my %depth_cols = map{ chomp; split( /\t/ ) }`cat $depth_col_file`;


$input_filename = 'data_mutations_extended.txt' if ( !defined $input_filename );


############# Find MAF files recursively if input is a directory ##################

if ( -d $input ){

    print "\nInput directory is $input\n\n";

    # Recursively find maf files for VEP annotation
    my @mafs = `find $input -name $input_filename`;
    foreach ( @mafs ){
        chomp;
        next if( !defined $_ );
        my ( $name, $path ) = fileparse( $_ );

        my $outFile;
        if ( !$output_maf ) {
            $outFile = ConstructOutputFileName( $_ );
        } else {
            $outFile = $path.$output_maf;
        }
        AnnotateMAF( $_, $outFile, $path.$annotated_maf );
    }
}
else{
    
    print "\nInput file is $input\n";

    my $outFile;
    if ( !$output_maf ) {
        $outFile = ConstructOutputFileName( $input );
    } else {
        $outFile = $output_maf;
    }

    print "\nOutput file is $outFile\n\n";
    
    AnnotateMAF( $input, $outFile, $annotated_maf );
}

############################ Run VEP wrapper script #################################

sub AnnotateMAF {
    my ( $inFile, $outFile, $annFile ) = @_;
    my ( $name, $path ) = fileparse( $inFile );
    
    PreProcess( $inFile, "$inFile.preprocessed.maf", "$inFile.msic" );
    my $preInFile = $inFile;
    $inFile = "$inFile.preprocessed.maf";
    
    my $cols_para = GetColsToRetain( $inFile );
    my $cp_hd_cmd = "if [ -e $outFile ]; then head -5 $inFile | egrep \"^#\" > $outFile.tmp;egrep -v \"^#\" $outFile >> $outFile.tmp;mv $outFile.tmp $outFile; fi;rm $inFile";

    # Contruct command line
    my $vep_cmd;
    $vep_cmd  = "$vep_wrap_script --ref-fasta $ref_fasta --vep-path $vep_path --vep-data $vep_data --vep-forks $vep_forks --input-maf $inFile --output-maf $outFile";
    $vep_cmd .= " --custom-enst $custom_enst_file" if (defined $custom_enst_file);
    $vep_cmd .= " --tmp-dir $tmp_dir" if ( $tmp_dir );
    $vep_cmd .= " $cols_para";
    
    if ( defined $use_cluster and (lc( $use_cluster ) eq 'yes' || lc( $use_cluster ) eq 'y' )){
        # Submit a job to run VEP annotation
        my $NR = `wc -l < $inFile`;
        chomp $NR;
    
        my $SP = "";
        $SP = " -sp 100 " if( $NR >= 100000 );
        
        # Submit jobs to annotate MAF files
        # print "bsub $SP -n $vep_forks -R \047span[hosts=1]\047 -eo $path/vep.log 'perl $vep_cmd;$cp_hd_cmd'\n";
        print "$preInFile : ";
        print `bsub $SP -n $vep_forks -R \047span[hosts=1]\047 -eo $path/vep.log \047perl $vep_cmd;$cp_hd_cmd\047`;
    }else{
        # Run VEP interatively
        # print "perl $vep_cmd\n";
        print `perl $vep_cmd`;
        
        # Replace MAF header
        system( $cp_hd_cmd ) == 0 or die "\nERROR: Failed to copy maf file header to annotated file!\nCommand: $cp_hd_cmd\n";
    }
}

######################### Construct output file name ##############################

sub ConstructOutputFileName {
    my $inFile = shift;
    my ( $name, $path ) = fileparse( $inFile );
    
    if ( $name =~ /extended/ ) {
        $name =~ s/extended/vep/;
    }elsif ( $name =~ /txt$/ ){
        $name =~ s/txt/vep.txt/;
    }elsif ( $name =~ /maf$/ ){
        $name =~ s/maf/vep.maf/;
    }else{
        $name .= 'vep.txt';
    }
    my $outFile = $path.$name;
    return $outFile;
}

############################ Data sanity checking ###################################

sub PreProcess {
    
    my ( $inFile, $outFile, $miscFile ) = @_;

    # Detect and replace DOS/MAC line ending symbol
    my $preInFile = `file $inFile`;
    if ( $preInFile =~ /CRLF line terminators/ ){
        `tr -d \$\047\\r\047 < $inFile > $outFile.tmpMaf`;
        $preInFile = $inFile;
        $inFile = "$outFile.tmpMaf";
    } elsif ( $preInFile =~ /CR line terminators/ ){
        `tr \047\\r\047 \047\\n\047 < $inFile > $outFile.tmpMaf`;
        $preInFile = $inFile;
        $inFile = "$outFile.tmpMaf";
    } else {
        $preInFile = $inFile;
    }
    
    # Read header line
    my @fields = split( /\t/, `head -5 $inFile | egrep -i "Hugo_Symbol|Chromosome|Tumor_Sample_Barcode"` );
    
    my ($idx, $chr_idx, $pos_idx, $ref_allele_idx, $tumor_allele_idx, $tumor_barcode_idx) = (0,-1,-1,-1,-1,-1);
    foreach ( @fields ){
        $chr_idx = $idx             if(lc($_) eq "chromosome");
        $pos_idx = $idx             if(lc($_) eq "start_position");
        $ref_allele_idx = $idx      if(lc($_) eq "reference_allele");
        $tumor_allele_idx = $idx    if(lc($_) eq "tumor_seq_allele2");
        $tumor_barcode_idx = $idx if(lc($_) eq "tumor_sample_barcode");
        $idx++;
    }
    
    if ($chr_idx==-1 or $pos_idx==-1 or $ref_allele_idx==-1 or $tumor_allele_idx==-1 or $tumor_barcode_idx==-1){
        die "Error: The MAF file lacks necessary column (see https://github.com/mskcc/vcf2maf/blob/master/data/minimalist_test_maf.tsv)\n";
    }
    
    open IN,  '<', $inFile  or die "Error: Unable to open $inFile\n";
    open OUT, '>', $outFile or die "Error: Unable to create $outFile\n";

    my $misc;
    while ( <IN> ) {
        next if ( /^\s*$/ );
        if (/^#|^Hugo_Symbol|^Chromosome|^Tumor_Sample_Barcode/i){
            print OUT $_;
        } else {
            my @F = split( /\t/ );
        
            if( defined $F[$chr_idx] and
                defined $F[$pos_idx] and
                defined $F[$ref_allele_idx] and
                defined $F[$tumor_allele_idx] and
                defined $F[$tumor_barcode_idx]){
                $F[$chr_idx] =~ s/23/X/;
                $F[$chr_idx] =~ s/24/Y/;
                print OUT join( "\t", @F );
            } else {
                $misc = $misc ? $misc . $_ : $_;
            }
        }
    }

    # Move illeagle mutations into $outFile.misc
    if (defined $misc) {
        open  MOUT, '>', $miscFile or die "Error: Unable to create $miscFile\n";
        print MOUT $misc;
        close MOUT;
    }

    close IN;
    close OUT;

    # Delete temperary file
    `rm $inFile` if ( $preInFile ne $inFile );
}

###################################################################################
# Construct a string for argumnet '--retain-cols' to keep columns in input MAF file

sub GetColsToRetain {
    my $myMAF = shift;

    my $ret = `head -5 $myMAF | egrep -i "Hugo_Symbol|Chromosome|Tumor_Sample_Barcode"`;
    chomp $ret;
    return '' if ( !defined $ret );

    my @columns = split( /\t/, $ret );
    
    my $cols2keep='--retain-cols';
    my $cnt=0;
    
    foreach ( @columns ){
        next if ( /ONCOTATOR/i );
        next if ( /Hugo_Symbol/i or /Entrez_Gene_Id/i or /Chromosome/i or /Start_Position/i or /End_Position/i or /NCBI_Build/i or /Variant_Classification/i or /Variant_Type/i );
        
        if ( exists $depth_cols{$_} ){
            $cols2keep = '--'.$depth_cols{ $_ }.' '.$_.' '.$cols2keep;
        }else{
            #s/ |\(|\)/_/g;
            if( $cnt==0 ){
                $cols2keep .= " \"$_";
                $cnt++;
            }else{
                $cols2keep .= ",$_";
            }
        }
    }
    return $cols2keep.'"';
}


__DATA__

=head1 NAME
 
 run-vep-wrapper - Annotate the effects of variants in MAF files.
 
=head1 SYNOPSIS
 
 perl run-vep-wrapper.pl --help
 perl run-vep-wrapper.pl -i /data/cbio-portal-data/studies/
 perl run-vep-wrapper.pl -i data_mutations_extended.txt -o data_mutations_vep.txt -t /tmp/wangq
 
=head1 OPTIONS
 
 -i | --input         A input MAF file or directory. If not specified, will look for it in a configuration file
 -o | --output-maf    Path to output MAF. If not specified, will use input filename appended with 'vep' as output filename
 -a | --annotated-maf Path to annotated MAF file. If not specified, will look for filename in a configuration file
 -c | --config-file   A configuration file. If not specified, will look for it in program, working, home directory, respectively
 -t | --tmp-dir       Folder for intermediate files. Override same variable in configuration file if specified
 -u | --use-cluster   Value: yes, no. 'yes' will submit job to cluster. 'no' runs vep interatively. Default is 'no'
 -v | --vep-forks     Number of CPUs to run VEP
 -h | -? | --help     Displays this information
 --man                Print detailed manual
 
=head1 DESCRIPTION
 
 The script provides application-level wrapper functions. The script firstly checks data sanity. When annotating a MAF, it keeps as many columns in the original MAF as possible. If a directory is provided as input, all maf files in the directory and its sub directories, whose names are specified in a configuration file, are recursively collected and annotated. The script can also submit jobs to computer cluster (using bsub command).
 
=head1 AUTHORS
 
 Qingguo Wang
 
 Nikolaus Schultz Lab
 Memorial Sloan Kettering Cancer Center
 New York, NY 10065
 
=head1 ACKNOWLEDGEMENTS
 
 Thank Cyriac Kandoth, Frederick Criscuolo, and Onur Sumer for suggestions and discussion
 
=cut