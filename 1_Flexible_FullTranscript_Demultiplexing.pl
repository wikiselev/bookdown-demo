use strict;
use warnings;

# Input a pair of sequencing FastQ files
# gunzip, read, write out smaller broken down files in a format suitable for submitting job array, and re-gzip each one in turn
# Breakdown by lane & cellID
# Keep order. 
# This should work equally well for single-end reads and can take any number of files as arguments.
# TESTED

if (@ARGV != 7) { 
	print STDERR "perl 1_Flexible_FullTranscript_Demultiplexing.pl read1.fq read2.fq b_pos b_len index mismatch prefix\n";
	print STDERR "
		read1.fq : barcode containing read
		read2.fq : non-barcode containg read
		b_pos : position of cell-barcode in the read. [\"start\" or \"end\"]
		b_len : length of cell-barcode (bp)
		index : file contain a single column of expected barcodes
		mismatch : maximum number of permitted mismatches (recommend 2)
		prefix : prefix for output fq files.\n";
	exit(1);
}

my $infile1 = $ARGV[0];
my $infile2 = $ARGV[1];
my $barcode_pos = $ARGV[2];
my $barcode_len = $ARGV[3];
my $barcode_index_file = $ARGV[4];
my $MAXmismatch = $ARGV[5];
my $OUTprefix = $ARGV[6];

if ($OUTprefix =~ /^(.+)\/[^\/]$/) {
	if ($1 ne ".") {
		system("mkdir -p $1");
	}
}

my %CellBarcodes = ();
my %ofhs1 = ();
my %ofhs2 = ();
open(my $ifh, $barcode_index_file) or die "Cannot open $barcode_index_file\n";
my $index=1;
while (<$ifh>) {
	chomp;
	$CellBarcodes{$_} = $index;
	open(my $fh1, '>', "$OUTprefix\_$_\_read1.fq") or die $!;
	$ofhs1{$index} = $fh1;
	open(my $fh2, '>', "$OUTprefix\_$_\_read2.fq") or die $!;
	$ofhs2{$index} = $fh2;
	$index++;
} close($ifh);


my $NotProperBarcodes = 0;
my $NotPossibleCell = 0;
my $AmbiguousCell = 0;
my $ExactMatch = 0;
my $Mismatch = 0;
my $total_reads = 0;
my $OutputReads = 0;

open(my $ifh1, $infile1) or die $!;
open(my $ifh2, $infile2) or die $!;
while(<$ifh1>) {
	my $file1line=$_;
	my $file2line = <$ifh2>;
	if ($file1line =~ /^@/) {
		# Ensure matching pair of reads
		my @thing1 = split(/\s+/, $file1line);
		my @thing2 = split(/\s+/, $file2line);
		my $readname = $thing1[0];
		#if ($readname ne $thing2[0]) {die "file1 & file2 readnames don't match!\n";}
		my $barcode_read = <$ifh1>;
		chomp $barcode_read;
		my $read2 = <$ifh2>;
		chomp $read2;
		$total_reads++;

		<$ifh1>; <$ifh2>;
		my $file1qual = <$ifh1>;
		chomp $file1qual;
		my $file2qual = <$ifh2>;
		chomp $file2qual;
		my $CellID = "";
		if ($barcode_pos eq "start") {
			$CellID = substr($barcode_read, 0, $barcode_len, "");
			substr($file1qual, 0, $barcode_len, "");
		} else {
			$CellID = substr($barcode_read, -$barcode_len, $barcode_len, "");
			substr($file1qual, -$barcode_len, $barcode_len, "");
		}
		my $mismatches = 0;
		if (!exists($CellBarcodes{$CellID})) {
			my @matches = ();
			my %close = ();
			foreach my $expected_barcode (keys(%CellBarcodes)) {
				my $count = ( $expected_barcode ^ $CellID ) =~ tr/\0//;
				if ($count >= length($expected_barcode)-$MAXmismatch) {
					$close{$expected_barcode} = $count;
				}
			}
			if (scalar(keys(%close)) > 0) {
				my $max = my_max(values(%close));
				$mismatches = length($CellID) - $max;
				foreach my $code (keys(%close)) {
					if ($close{$code} == $max) {
						push(@matches, $code);
					}
				}
			}
			if (scalar(@matches) == 1) {
				$CellID = $matches[0];
				$Mismatch++;
			} elsif (scalar(@matches) > 1) {
				$AmbiguousCell++;
				next;
			} else {
				$NotPossibleCell++;
				next;
			}
		} else {
			$ExactMatch++;
		}
		# print the read
		my $handle1 = $ofhs1{$CellBarcodes{$CellID}};
		my $handle2 = $ofhs2{$CellBarcodes{$CellID}};
		print $handle1 "$readname\n$barcode_read\n+\n$file1qual\n";
		print $handle2 "$readname\n$read2\n+\n$file2qual\n";
		$OutputReads++;
	} else {next;}
}

print STDERR "
Doesn't match any cell: $NotPossibleCell
Ambiguous: $AmbiguousCell
Exact Matches: $ExactMatch
Contain Mismatches: $Mismatch
Input Reads: $total_reads
Output Reads: $OutputReads\n";
close($ifh1);
close($ifh2);
foreach my $ofh1 (keys(%ofhs1)) {close($ofhs1{$ofh1});}
foreach my $ofh2 (keys(%ofhs2)) {close($ofhs2{$ofh2});}



sub my_max {
	if (scalar(@_) == 1) {return($_[0])};
	my $max = shift;
	foreach my $ele (@_) {
		if ($ele > $max) {$max = $ele;}
	}
	return($max);
}
