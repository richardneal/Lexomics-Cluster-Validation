#!/usr/bin/perl

use strict;
use warnings;

open DATATABLE, "AlwaysRightData.tsv" or die $!;
open SAVE, "AlwaysRightChunksShortened.txt" or die$!;
open OUTPUT, ">", "AlwaysRightData5and5.tsv" or die$!;

my @toSave = <SAVE>;
chomp(@toSave);

my $firstLine = <DATATABLE>;
print OUTPUT ($firstLine);

my @labelUsed;

for my $i (0 .. @toSave - 1)
{
   $labelUsed[$i] = 0;
}

while (<DATATABLE>) {
	my $line = $_;
	$line =~ /([\w|-]*)\t.*/;
        #print($1);

        for my $i (0 .. @toSave - 1)
        {
        	#print $_;
        	if(index($1,$toSave[$i]) != -1)
                {
                        $labelUsed[$i] = 1;
                        #print $_;
                	print OUTPUT ($line);
                }
        }

}

print("Done");

#for my $i (0 .. @toSave - 1)
#{
#   if(!$labelUsed[$i])
#   {
#   	print($i);
#	print("\n");
#   }
#}

close DATATABLE ;
close SAVE;
close OUTPUT;