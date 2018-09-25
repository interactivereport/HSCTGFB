#!/usr/bin/perl

# Author: Baohong Zhang (baohong.zhang@pfizer.com)
# It is used to filter out very low count, if more than 80% of samples has read count less than $MINREAD and the average count on all samples
# not greater than user specified number (e.g. 2), then the gene will be filtered out.

use strict;

if (@ARGV < 3) { 
	print "Usage: $0 <count file> <min avg count> <count starting column> <file of exclude genes>\n";
	print "e.g. $0 IL6-Lupus_RNA-seq_count_479-US-1483814.txt 2 3 ~/pipeline/RNA-seq/globin_exclude.txt > IL6-Lupus_RNA-seq_count_479-US-1483814_flt.txt\n";

	exit;
}



my $MINAVG = $ARGV[1];

my ($start, $i, $line, @items, $n, $reads, $tiny);

my $MINREAD = 10;
my $MINREAD_PERCENT = 0.8;  # one of the condition for filtering: genes with less than $MINREAD for more than 80% samples.


my %geneExclude = ();
if (open (FILE, "$ARGV[3]")) {
	while($line = <FILE>) {
	    next if ($line =~ m/^#/);
	    $line =~ s/\s+$//g;
	    $geneExclude{$line} = 1;
	}
	close(FILE);
}

open (FILE, "$ARGV[0]") || die $@;

$line = <FILE>;
$line =~ s/\s+$//g;
@items = split(/\t/, $line);

print $items[0];
for ($i=1; $i<@items; $i++) {
	print "\t$items[$i]";
}
print "\n";

$start = $ARGV[2] - 1;

while ($line = <FILE>) {
	$line =~ s/\s+$//g;
	@items = split(/\t/, $line);

	$n = 0;
	$tiny = 0;
	$reads = 0;
	for ($i=$start; $i<@items; $i++) {
		$n += 1;
		$reads += $items[$i];

		if ($items[$i] < $MINREAD) {
			$tiny += 1;
		}
	}

	if (($tiny / $n) > $MINREAD_PERCENT && ($reads/$n) <= $MINAVG) {
		next;
	}

	next if (defined $geneExclude{$items[0]});

	$items[0] =~ s/\.\d+$//g;
	print $items[0];
	for ($i=1; $i<@items; $i++) {
		print "\t$items[$i]";
	}
	print "\n";
}
close(FILE);
