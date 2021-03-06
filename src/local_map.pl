#!/usr/bin/perl
######################################################################
#                            local_map                               #
######################################################################
#
#     cegma pipeline code
#
######################################################################
#
# 	$Id: local_map.pl,v 1.5 2011/07/19 17:31:12 keith Exp $

use strict; use warnings; use sigtrap;
use FAlite;
use Getopt::Std;
use vars qw($opt_p $opt_h $opt_n $opt_g $opt_d $opt_f $opt_i $opt_v);

##########################################################################
##                                                                      ##
##                     VARIABLES  AND DEFINITIONS                       ##
##                                                                      ##
##########################################################################

my $CEGMA = defined( $ENV{'CEGMA'} )    ? $ENV{'CEGMA'}    : '.';
my $TMP   = defined( $ENV{'CEGMATMP'} ) ? $ENV{'CEGMATMP'} : '/tmp';

getopts('p:n:h:g:d:fi:v');

my $PREFIX = "";
my $HMMDIR;
my $PROTOPT;
my $WISE_GFF;
my $GENEID     = "$CEGMA/data/generic.param";
my $HMM_PREFIX = "KOG";

$PREFIX     = $opt_n if $opt_n;
$HMMDIR     = $opt_h if $opt_h;
$PROTOPT    = $opt_p if $opt_p;
$WISE_GFF   = $opt_g if $opt_g;
$GENEID     = $opt_d if $opt_d;
$HMM_PREFIX = $opt_i if $opt_i;

die "

# This script takes a protein multifasta file and a genomic multifasta file
# and runs genewise and geneid+genewise and evaluate the results.

# we considered the sequence of the proteins and the genomic to be 
# in the same order.

local_map.pl <genomic_fasta_file>
options:
   -n <prefix of the output files >
   -p <prot_fasta_file>
   -h <directory containing the hmm profiles>
   -g <existing genewise_gff file>    # This option skips genewise  
   -d <geneid_parameter>
   -i <prefix to use for hmm output
   -f force the prediction of one gene structure with geneid
   -v verbose mode, show ID of each KOG that is being processed

" unless ( @ARGV == 1 && ( $opt_h || $opt_p ) );

# Bash variable:
my $genome_file = $ARGV[0];

# Security check to not overwrite the results
die "$genome_file does not exist\n" if ( not -e $genome_file );
die "Can't find $WISE_GFF\n" if ( defined $WISE_GFF and not -e $WISE_GFF );
print "$PREFIX.genewise exists and will be overwritten\n" if ( -e "$PREFIX.genewise" and ( not defined $WISE_GFF ) );
print "NOTE: Will use specifed $WISE_GFF file instead of running genewise\n" if ( defined $WISE_GFF );

# remove existing files else they will just be added to rather than replaced
if ( -e "$PREFIX.geneid" ) {
    print "WARNING: removing existing version of $PREFIX.geneid\n";
    unlink("$PREFIX.geneid") or die "Can't remove $PREFIX.geneid\n";
}
if ( -e "$PREFIX.geneid.gff" ) {
    print "WARNING: removing existing version of $PREFIX.geneid.gff\n";
    unlink("$PREFIX.geneid.gff") or die "Can't remove $PREFIX.geneid.gff\n";
}

# Proteins option
my $PROT;

if ($PROTOPT) {
    open( PROT_FILE, "<", $PROTOPT ) or die "Can't open $PROTOPT\n";
    $PROT = new FAlite( \*PROT_FILE );
}

# Opening the fasta files
open( GENOMIC_FILE, "<", $genome_file ) or die "Can't open $genome_file\n";
my $GENOMIC = new FAlite( \*GENOMIC_FILE );

# Main loop of execution
while ( my $gen_entry = $GENOMIC->nextEntry ) {

    # Opening temporary fasta file, which will contain sequence of one 'chunk'
    open( TMP_GENOMIC, ">", "$TMP/genomic$$.fa" ) or die "Can't write to $TMP/genomic$$.fa\n";
    print_fasta( $gen_entry, \*TMP_GENOMIC );
    my ($kog_id)     = $gen_entry->def =~ /($HMM_PREFIX[0-9]+)/;
    my ($kog_id_sub) = $gen_entry->def =~ /($HMM_PREFIX[0-9.]+)/;
    close(TMP_GENOMIC);

    print "Processing chunk: $kog_id_sub\n" if ($opt_v);

# if we are using a specified protein file instead of using a HMM for genewise, need to write that to a file
    if ($PROTOPT) {
        open( TMP_PROT, ">", "$TMP/prot$$.fa" ) or die "Can't write to $TMP/prot$$.fa\n";
        my ($prot_entry) = $PROT->nextEntry;
        print_fasta( $prot_entry, \*TMP_PROT );
        close(TMP_PROT);
    }

    # do we need to run Genewise?
    if ( !$WISE_GFF ) {

        # form command based on whether $PROTOPT has been specified
        my $command = "genewise -splice_gtag -quiet -gff -pretty -alb -hmmer $HMMDIR/$kog_id.hmm $TMP/genomic$$.fa > $TMP/genewise$$";
        $command = "genewise -splice_gtag -quiet -gff -pretty -alb $TMP/prot$$.fa $TMP/genomic$$.fa > $TMP/genewise$$" if ($PROTOPT);

        system("$command") && die "Can't run $command\n";

        # convert genewise file to GFF
        system("$CEGMA/bin/parsewise $TMP/genewise$$ > $TMP/genewise$$.gff") && die "Can't run parsewise\n";

        # add tmp results to main genewise output files
        system("cat $TMP/genewise$$     >> $PREFIX.genewise") && die "Can't concatenate $TMP/genewise$$ to $PREFIX.genewise\n";
        system("cat $TMP/genewise$$.gff >> $PREFIX.genewise.gff") && die "Can't concatenate $TMP/genewise$$.gff to $PREFIX.genewise.gff\n";

   # not sure why we also add this to the main genewise GFF file for each entry
   #		system ("echo \"#\$\" >> $PREFIX.genewise.gff") && die "Can't run echo\n";
    } else {

        # use pre-existing genewise file provided
        system("grep -w $kog_id_sub $WISE_GFF > $TMP/genewise$$.gff");
    }

# at this point the $TMP/genewise$$.gff file might be empty if there were no genewise results
# maybe we could skip the following steps if this is the case?
    unless ( -z "$TMP/genewise$$.gff" ) {

        # Obtanining the boundaries of the genewise alignment
        open( TMP_GENEWISE, "<", "$TMP/genewise$$.gff" ) or die "Can't open $TMP/genewise$$.gff\n";

        my $initial_genewise_coord = 0;
        my $final_genewise_coord   = 0;
        my $counter                = 0;

        while (<TMP_GENEWISE>) {
            my @gff = split;
            $initial_genewise_coord = $gff[3] if ( $counter == 0 );
            $counter++;
            $final_genewise_coord = $gff[4];
        }
        close(TMP_GENEWISE);

# now want to extract just the region of the genewise match from the genome (+/- 1000 nt)
# and also mask out some regions with N characters

        open( TMP_GENOMIC, ">", "$TMP/genomic_mask$$.fa" ) or die "Can't write to $TMP/genomic_mask$$.fa\n";
        print_fasta_masked( $gen_entry, \*TMP_GENOMIC, $initial_genewise_coord - 1000, $final_genewise_coord + 1000 );
        close(TMP_GENOMIC);

        # running geneid
        my $gff_command = "geneid  -P $GENEID -S $TMP/genewise$$.gff -WG $TMP/genomic_mask$$.fa | sed 's/exon/Exon/' >> $PREFIX.geneid.gff";
        my $fa_command = "geneid -AP $GENEID -S $TMP/genewise$$.gff -W  $TMP/genomic_mask$$.fa | sed 's/exon/Exon/' >> $PREFIX.geneid";

        # use force single prediction mode? (-F)
        if ($opt_f) {
            $gff_command = "geneid  -P $GENEID -S $TMP/genewise$$.gff -WGF $TMP/genomic_mask$$.fa | sed 's/exon/Exon/' >> $PREFIX.geneid.gff";
            $fa_command = "geneid -AP $GENEID -S $TMP/genewise$$.gff -WF  $TMP/genomic_mask$$.fa | sed 's/exon/Exon/' >> $PREFIX.geneid";
        }
        system("$gff_command") && die "Can't run $gff_command\n";
        system("$fa_command")  && die "Can't run $fa_command\n";
    }

   # do we need these?
   #    system ("echo \"#\$\" >> $PREFIX.geneid.gff") && die "Can't run echo\n";
   #    system ("echo \"#\$\" >> $PREFIX.geneid") && die "Can't run echo\n";
}

close(PROT_FILE) if ($PROTOPT);
close(GENOMIC_FILE);

system("awk '{if (substr(\$1,1,1)==\">\") bol=1;if (\$0==\"\") bol=0;if (bol) print; }' $PREFIX.geneid > $PREFIX.geneid.fa") && die "Can't run awk\n";

# removing temporary files
my @files = ( "$TMP/genomic$$.fa", "$TMP/genomic_mask$$.fa", "$TMP/genewise$$", "$TMP/genewise$$.gff" );
push( @files, "$TMP/prot$$.fa" ) if ($PROTOPT);

foreach my $file (@files) {
    unlink("$file") or warn "WARNING: Can't remove $file\n";
}

exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################
##########################################################################

##########################################################################
##                                                                      ##
##                            SUBROUTINES                               ##
##                                                                      ##
##########################################################################

# simple subroutine to print a FASTA file
sub print_fasta {

    my ( $entry, $file ) = @_;
    my ( $id,    $seq )  = split( /\n/, $entry );
    my ( $i,     $e )    = ( 1, length($seq) );
    print $file "$id\n";
    while ( $i <= $e ) {
        print $file substr( $seq, $i - 1, 60 ) . "\n";
        $i += 60;
    }
}

# This subrutine masks with N from 0 to pos1 and from pos2 to the end
sub print_fasta_masked {

    my ( $entry, $file, $pos1, $pos2 ) = @_;

    my ( $id, $seq ) = split( /\n/, $entry );
    my ( $i, $e ) = ( 1, length($seq) );

    ( $pos1 = 0 )  if ( $pos1 < 0 );
    ( $pos2 = $e ) if ( $pos2 <= $pos1 );

    # printer header
    print $file "$id\n";

  # this doesn't work properly if there is no genewise hit,
  # $pos 2 ends up = 1000 which means first 1000 nt are not masked
  #	print "i = $i, e = $e, pos1 = $pos1, pos2 = $pos2, L = ", length($seq),"\n";

    while ( $i < $pos1 ) {
        print $file "N";
        print $file "\n" if ( !( $i % 60 ) );
        $i++;
    }

    while ( $i <= $e && $i <= $pos2 ) {
        print $file substr( $seq, $i - 1, 1 );
        print $file "\n" if ( !( $i % 60 ) );
        $i++;
    }

    while ( $i < $e ) {
        print $file "N";
        print $file "\n" if ( !( $i % 60 ) );
        $i++;
    }
}
