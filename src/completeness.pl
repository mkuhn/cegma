#!/usr/bin/perl
######################################################################
#                           Completeness                             #
######################################################################
#
#     cegma pipeline code
#
######################################################################

# 	$Id: completeness.pl,v 1.5 2011/07/19 17:31:12 keith Exp $

use strict; use warnings; use sigtrap;
use FAlite;
use Getopt::Std;
use vars qw($opt_m $opt_e $opt_v);

getopts('mev');

die "

# This script takes a hmmsearch output file and a cutoff paramaters file 
# and give the statistic of completeness

hmm_select.pl <hmmsearch_output> <selected KOG ids> <cutoff_file>

options:

  -m           printing missing genes
  -v           verbose

" unless ( @ARGV == 3 );

# want to process just the one top scoring gene prediction in each KOG
my $hmm_file           = $ARGV[0];
my $selected_kogs_file = $ARGV[1];
my $cutoff_file        = $ARGV[2];

######################################################################
#
#     Security check
#
######################################################################

if ( !( -e "$hmm_file" ) ) {
    print "File does not exist: $hmm_file \n";
    exit(1);
}

if ( !( -e "$selected_kogs_file" ) ) {
    print "File does not exist: $selected_kogs_file \n";
    exit(1);
}

if ( !( -e "$cutoff_file" ) ) {
    print "File does not exist: $cutoff_file \n";
    exit(1);
}

######################################################################
#
#      Reading the cutoff and conservation group for each protein
#
######################################################################

my %score;
my %conservation;
my %conservation_groups;
my %score_cutoff;
my %prot_length;
my $total_proteins;

my %length_cutoff;

if ($opt_e) {
    %length_cutoff = (
        "Cutoff 90%" => '90',
        "Cutoff 80%" => '80',
        "Cutoff 70%" => '70',
        "Cutoff 60%" => '60',
        "Cutoff 50%" => '50',
        "Cutoff 0%"  => '0',
    );

} else {
    %length_cutoff = (
        Complete => '70',
        Partial  => '0',
    );
}

open( CUTOFF_FILE, "<", "$cutoff_file" ) or die "Can't read from $cutoff_file\n";
while (<CUTOFF_FILE>) {
    if (/(.*) (.*) (.*) (.*)/) {
        $conservation{$1} = $2;
        $conservation_groups{$2}++;
        $score_cutoff{$1} = $3;
        $prot_length{$1}  = $4;
        $score{$1}        = 0;
        $total_proteins++;
    }
}

close(CUTOFF_FILE);

######################################################################
#
#     Reading the selected KOGs id file to get IDS of just the best
#     gene prediction in each KOG
#
######################################################################

my %best_kog_ids;

open( IDS, "<", "$selected_kogs_file" ) or die "Can't read from $selected_kogs_file\n";
while ( my $id = <IDS> ) {
    chomp($id);
    $best_kog_ids{$id} = 1;
}
close(IDS);

######################################################################
#
#     Reading the hmmsearch output and getting the score and length
#     of the alignments
#
######################################################################

my %over_cutoff;
my %over_cutoff_locus;
my %over_cutoff_group;
my %over_cutoff_group_locus;
my $locus;
my $warning_version = 0;

open( HMMSEARCH, "<", "$hmm_file" ) or die "Can't read from $hmm_file\n";

# check HMMER version
<HMMSEARCH>;
my $second_line = <HMMSEARCH>;
my ($version) = $second_line =~ m/HMMER (\d+\.\d+)/;
$warning_version = 1 if ( $version ne "3.0" );

while (<HMMSEARCH>) {

    next if m/^#/;
    if (m/^>>/) {

        # grab specific KOG ID
        #		>> KOG0328.15_1|geneid_v1.4_predicted_protein_1|408_AA
        my ($kog_id) = $_ =~ m/>>\s+(\S+?)_/;
        my $id = $kog_id;
        $id =~ s/\..*//;
        next unless exists( $best_kog_ids{$kog_id} );

        print "Evaluating $kog_id" if ($opt_v);

        if ( defined $score_cutoff{$id} ) {
            print "\tPART OF 248 CEGs set\n" if ($opt_v);
        } else {
            print "\n" if ($opt_v);
            next;
        }

        # skip 2 more lines
        <HMMSEARCH>; <HMMSEARCH>;
        my $score_line = <HMMSEARCH>;
        my ( $index, undef, $score, undef, undef, undef, $start, $end ) = $score_line =~ m/\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;

        if ( $index == 1 && $score > $score_cutoff{$id} ) {
            my $align_length = $end - $start + 1;
            for my $cutoff ( keys %length_cutoff ) {
                my $target = ( $prot_length{$id} / 100 * $length_cutoff{$cutoff} );

                if ( $align_length > ( $prot_length{$id} / 100 * $length_cutoff{$cutoff} ) ) {
                    $over_cutoff{$cutoff}++;
                    $over_cutoff_locus{$cutoff}{$id}++;
                    $over_cutoff_group{$cutoff}{ $conservation{$id} }++;
                    $over_cutoff_group_locus{$cutoff}{ $conservation{$id} }{$id}++;
                }
            }
        }
    }
}

close(HMMSEARCH);

######################################################################
#
#     Computing averages and other statistics
#
######################################################################

my %average;
my %average_group;
my %locus_orth;
my %locus_orth_group;
my %percent_orth;
my %percent_orth_group;

for my $cutoff ( sort keys %length_cutoff ) {

    if ( $over_cutoff{$cutoff} ) {
        $average{$cutoff} = $over_cutoff{$cutoff} / ( scalar keys %{ $over_cutoff_locus{$cutoff} } );
    } else {
        $average{$cutoff} = 0;
    }

    for my $locus ( keys %{ $over_cutoff_locus{$cutoff} } ) {
        if ( $over_cutoff_locus{$cutoff}{$locus} > 1 ) {
            $locus_orth{$cutoff}++;
        }
    }

    if ( $over_cutoff{$cutoff} && $locus_orth{$cutoff} ) {
        $percent_orth{$cutoff} = $locus_orth{$cutoff} / ( scalar keys %{ $over_cutoff_locus{$cutoff} } ) * 100;
    } else {
        $percent_orth{$cutoff} = 0;
    }

    for my $group ( sort keys %conservation_groups ) {

        if ( $over_cutoff_group{$cutoff}{$group} ) {
            $average_group{$cutoff}{$group} = $over_cutoff_group{$cutoff}{$group} / ( scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} } );
        } else {
            $average_group{$cutoff}{$group} = 0;
        }

        for my $locus ( keys %{ $over_cutoff_group_locus{$cutoff}{$group} } ) {
            if ( $over_cutoff_group_locus{$cutoff}{$group}{$locus} > 1 ) {
                $locus_orth_group{$cutoff}{$group}++;
            }
        }

        if ( $over_cutoff_group{$cutoff}{$group} && $locus_orth_group{$cutoff} && ( scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} } ) ) {
            $percent_orth_group{$cutoff}{$group} = $locus_orth_group{$cutoff}{$group} / ( scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} } ) * 100;
        } else {
            $percent_orth_group{$cutoff}{$group} = 0;
        }

    }

}

######################################################################
#
#     Printing the output
#
######################################################################

if ($warning_version) {
    print "\n WARNING executed HMMER version has not been tested !!! \n"
}

print "\n";

print "#      Statistics of the completeness of the genome based on $total_proteins CEGs      #\n\n";

print "              #Prots  %Completeness  -  #Total  Average  %Ortho \n";

print "\n";

for my $cutoff ( sort keys %length_cutoff ) {

    if ( scalar keys %{ $over_cutoff_group_locus{$cutoff} } && $over_cutoff{$cutoff} ) {
        printf "%10s %8d %11.2f      - %5d %8.2f %9.2f\n",
          $cutoff,
          scalar keys %{ $over_cutoff_locus{$cutoff} },
          ( scalar keys %{ $over_cutoff_locus{$cutoff} } ) / $total_proteins * 100,
          $over_cutoff{$cutoff},
          $average{$cutoff},
          $percent_orth{$cutoff};
    }
    else {
        printf "%10s %8d %11.2f      - %5d %8.2f %9.2f\n",
          $cutoff, 0, 0, 0, 0, 0, 0;

    }

    print "\n";

    for my $group ( sort keys %conservation_groups ) {
        if ( scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} } ) {
            printf "   Group %s %8d %11.2f      - %5d %8.2f %9.2f\n",
              $group,
              scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} },
              ( scalar keys %{ $over_cutoff_group_locus{$cutoff}{$group} } ) / $conservation_groups{$group} * 100,
              $over_cutoff_group{$cutoff}{$group},
              $average_group{$cutoff}{$group},
              $percent_orth_group{$cutoff}{$group};
        }
        else {
            printf "   Group %s %8d %11.2f      - %5d %8.2f %9.2f\n",
              $group, 0, 0, 0, 0, 0, 0;
        }
    }
    print "\n";
}

print "#    These results are based on the set of genes selected by Genis Parra   #\n\n";

if ($opt_m) {
    print "##              Missing proteins            ##\n";
    for my $cutoff ( sort keys %length_cutoff ) {
        print "# $cutoff \n";

        for my $locus ( sort keys %score_cutoff ) {
            if ( !( defined $over_cutoff_locus{$cutoff}{$locus} ) ) {
                print "$locus\n";
            }
        }
    }
}

exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################
