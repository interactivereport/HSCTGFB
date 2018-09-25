#!/usr/bin/perl

# Author: Baohong Zhang (baohong.zhang@pfizer.com)
# It is used to filter out very low value, if more than 80% of samples has RPKM value less than $MINVALUE and the average value on all samples
# not greater than user specified number (e.g. 2), then the gene will be filtered out.

use strict;
my ($start, $i, $line, @items, $n, $values, $tiny);
my $MINVALUE = 1;
my $MINVALUE_PERCENT = 0.8;  # one of the condition for filtering: genes with less than $MINVALUE for more than 80% samples.
my %idx = ();
my %sampleTissue = ();
my @colname = ();

if (@ARGV < 4) { 
	print "Usage: $0 <value file> <min avg value> <value starting column> <tissue> <file of exclude genes>\n";
	print "e.g. $0 IL6-Lupus_RNA-seq_value_479-US-1483814.txt 2 3 ~/pipeline/RNA-seq/globin_exclude.txt > IL6-Lupus_RNA-seq_value_479-US-1483814_flt.txt\n";

	exit;
}
my $MINAVG = $ARGV[1];
my $TISSUE = $ARGV[3];

open (FILE, "sample.annotation.txt") || die $@;
$line = <FILE>;
chop($line);
@items = split(/\t/, $line);
for ($i=0; $i<@items; $i++) {
	$idx{$items[$i]} = $i;
}

while ($line = <FILE>)  {
	$line = <FILE>;
	chop($line);
	@items = split(/\t/, $line);

	$sampleTissue{$items[$idx{'sample_id'}]} = $items[$idx{'Tissue'}];
}
close(FILE);




my %geneExclude = ();
if (open (FILE, "$ARGV[4]")) {
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
$colname[0] = $items[0];
for ($i=1; $i<@items; $i++) {
	$colname[$i] = $items[$i];	
	print "\t$items[$i]";
}
print "\n";

$start = $ARGV[2] - 1;

while ($line = <FILE>) {
	$line =~ s/\s+$//g;
	@items = split(/\t/, $line);

	$n = 0;
	$tiny = 0;
	$values = 0;
	for ($i=$start; $i<@items; $i++) {
		next if ($sampleTissue{$colname[$i]} =~ m/$TISSUE/i);
		$n += 1;
		$values += $items[$i];

		if ($items[$i] < $MINVALUE) {
			$tiny += 1;
		}
	}

	# lowly expressed genes are filtered out
print "\n\n===", $tiny / $n, "\t", $values/$n, "===\n";
	if (($tiny / $n) > $MINVALUE_PERCENT && ($values/$n) <= $MINAVG) {
		next;
	}

	next if (defined $geneExclude{$items[0]});

#	$items[0] =~ s/\.\d+$//g;
	print $items[0];
	for ($i=1; $i<@items; $i++) {
		if ($i >= $start) {
			next if ($sampleTissue{$colname[$i]} =~ m/$TISSUE/i);
		}
		print "\t$items[$i]";
	}
	print "\n";
}
close(FILE);
