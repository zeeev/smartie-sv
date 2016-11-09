#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

annotate_size_dist.pl -f my.fai -s my.svs.bed

Description:

Adds two additional annotations to a printGaps file:
1. The size of the contig.
2. The the distance of the SV to the beggining or end.
   Whichever is smaller.

";


my ($help);
my $fai;
my $svs;
my $opt_success = GetOptions('help'      => \$help,
                             "file=s"    => \$fai,
                             "svs=s"     => \$svs,
    );

die $usage if $help || ! $opt_success;

die $usage unless defined $fai && defined $svs;
open (my $IN, '<', $fai) or die "Can't open $fai for reading\n$!\n";

my %query_size;

while (<$IN>) {
    chomp;
    my @line = split /\t/, $_;
    $query_size{$line[0]} = $line[1];
}

close $IN;

open (my $INB, '<', $svs) or die "Can't open $svs for reading\n$!\n";

while(<$INB>){
    chomp;
    my @line = split /\t/, $_;
    my $seqid = $line[-3];
    my $startD = $line[-2] - 0;
    my $endD   = $query_size{$seqid} -  $line[-1];

    my $d = 0;

    $startD > $endD ? $d = $endD : $d = $startD;

    my $per = $d / $query_size{$seqid};

    push @line, $query_size{$seqid}, $d, $per;
    print join "\t", @line, "\n";
}

close $INB;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

