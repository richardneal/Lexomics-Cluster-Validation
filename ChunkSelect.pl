#!/usr/bin/perl

use strict;
use warnings;

#print("yes");

open DATATABLE, "noDup_10K_cull_merged_4mer.tsv" or die $!;
open SAVE, "OrganismsToSave.txt" or die$!;
open OUTPUT, ">", "NewGenomicsData.tsv" or die$!;

my @toSave = <SAVE>;
chomp(@toSave);

my $firstLine = <DATATABLE>;
print OUTPUT ($firstLine);

while (<DATATABLE>) {
	my $line = $_;
	$line =~ /(\w*)\t.*/;
        #print($1);

        foreach(@toSave)
        {
        	#print $_;
        	if(index($1,$_) != -1)
                {
                        print $_;
                	print OUTPUT ($line);
                }
        }
}

close DATATABLE ;
close SAVE;
close OUTPUT;