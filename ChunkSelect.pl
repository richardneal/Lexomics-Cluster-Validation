#!/usr/bin/perl

use strict;
use warnings;

open DATATABLE, "noDup_10K_cull_merged_4mer.tsv" or die $!;
open SAVE, "OrganismsToSave.txt" or die$!;
open OUTPUT, ">", "NewGenomicsData.tsv" or die$!;

my @toSave = <SAVE>;

my $firstLine = <DATATABLE>;
print OUTPUT ($firstLine);

while (<DATATABLE>) {
	my $line = $_;
	$line =~ /(\w*)\t.*/;
        #print($1);

        foreach(@toSave)
        {
        	if(index($1,$_) != -1)
                {
                	print OUTPUT ($line);
                }
        }
}

close DATATABLE ;
close SAVE;
close OUTPUT;