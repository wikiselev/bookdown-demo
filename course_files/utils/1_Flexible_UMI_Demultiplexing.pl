use strict;
use warnings;

# Matches upto two mismatches between observed cell barcodes and the expected cell barcodes
# Excludes reads with problematic UMIs: >= 80% A, >= 80% T, contained in adaptor sequence. - Note do not provide adaptors for short UMI datasets (fewer than 7 bases) since there is a high probability of real UMIs being contained in the adaptor for such cases.
# Allows barcodes to contain ambiguous bases
# Allows trailing bases at the end of the barcode sequence but requires barcodes to begin from the first base in the barcode sequence.

if (@ARGV != 6) { 
print STDERR "perl 1_Flexible_UMI_Demultiplexing.pl read1.fq read2.fq b_structure index mismatch prefix\n";
print STDERR "
		read1.fq : barcode/umi containing read
		read2.fq : non-barcode containing read
		b_structure : a single string of the format C##U# or U#C## 
			where C## is the cell-barcode and U# is the UMI.
			e.g. C10U4 = a 10bp cell barcode followed by a 4bp UMI
		index : file containg a single column of expected cell-barcodes.
			if equal to \"UNKNOWN\" script will output read counts for each unique barcode.
		mismatch : maximum number of permitted mismatches (recommend 2)
		prefix : prefix for output fastq files.\n";
exit(1);}
my $infile1 = $ARGV[0];
my $infile2 = $ARGV[1];
my $barcodestructure = $ARGV[2];
my $MAXmismatch = $ARGV[4];

# Parse Barcode Structure #



my $order = -1;
my $C_len = -1;
my $U_len = -1;

if ($barcodestructure =~ /^C(\d+)U(\d+)$/) {
	$order=1;
	$C_len = $1;
	$U_len = $2;
	print "Barcode Structure: $C_len bp CellID followed by $U_len bp UMI\n";
} elsif ($barcodestructure =~ /^U(\d+)C(\d+)$/) {
	$order = 0;
	$C_len = $2;
	$U_len = $1;
	print "Barcode Structure: $U_len bp UMI followed by $C_len bp CellID\n";
} else {
	die "$barcodestructure not recognized.\n";
}
# ----------------------- #

my $OUTprefix = $ARGV[5]; #prefix for output
#Ensure output directory exists
if ($OUTprefix =~ /^(.+)\/[^\/]$/) {
	if ($1 ne ".") {
		system("mkdir -p $1");
	}
}

# Read Expected Cell Barcodes #
my %CellBarcodes = ();
my %ofhs = ();
if ($ARGV[3] ne "UNKNOWN") {
	open (my $ifh, $ARGV[3]) or die "Cannot open $ARGV[3]\n";
	my $index=1;
	while (<$ifh>) {
		chomp;
		if ($_ =~/^#/) {next;}
		my $barcode = $_;
		$CellBarcodes{$barcode} = $index;
		open(my $fh,'>',"$OUTprefix\_$barcode.fq") or die $!;
		$ofhs{$index} = $fh;
		$index++;
	} close($ifh);
}
# --------------------------- #


### Process Reads ###

# Summary Statistics
my $NotProperBarcodes = 0;
my $NotPossibleCell = 0;
my $AmbiguousCell = 0;
my $ExactMatch = 0;
my $Mismatch = 0;
my $BadUMI = 0;
my $total_reads = 0;
my $OutputReads=0;

open (my $ifh1, $infile1) or die $!;
open (my $ifh2, $infile2) or die $!;
while(<$ifh1>) {
	my $file1line = $_;
	my $file2line = <$ifh2>;
	if ($file1line =~ /^@/) { #Skip any file headers

		# Ensure matching pair of reads
		my @thing1 = split(/\s+/,$file1line);
		my @thing2 = split(/\s+/,$file2line);
		#my $readname1 = chop($thing1[0]);
		#my $readname2 = chop($thing2[0]);
		#if ($readname1 ne $readname2) {die "file1 & file2 readnames don't match! $thing1[0] $thing2[0]\n";}
		my $readname = $thing1[0];
		my $barcodes = <$ifh1>;
		my $read = <$ifh2>;
		$total_reads++;

		# Parse barcodes
		my $CellID = ""; my $UMI = "";
		if ($order) {
			if ($barcodes =~ /^([ATCGNUKMRYSWBVHDX]{$C_len})([ATCGNUKMRYSWBVHDX]{$U_len})/) {
				$CellID = $1; $UMI = $2;
			} else {$NotProperBarcodes++; next;}
		} else {
			if ($barcodes =~ /^([ATCGNUKMRYSWBVHDX]{$U_len})([ATCGNUKMRYSWBVHDX]{$C_len})/) {
				$CellID = $2; $UMI = $1;

			} else {$NotProperBarcodes++; next;}
		}
#		if ($barcodes =~ /^([ATCGNUKMRYSWBVHDX]{11})([ATCGNUKMRYSWBVHDX]{10})/) {


		# Correct for upto two mismatches between observed and expected cell barcodes
		if ($ARGV[3] ne "UNKNOWN") {
		my $mismatches = 0;
		if (!exists($CellBarcodes{$CellID})) { # Not an expected barcode

			# Barcode contains uncertain bases -> convert to wildcards and pattern match on expected barcodes. (given priority over barcodes with higher confidence mismatches)
			if ($CellID !~ /^[ATCG]+$/) {
				$mismatches = () = $CellID =~ /[^ATCG]/g; # count uncertain bases as mismatches
				$CellID =~ s/[^ATCG]/./g; #Turn non-ATCG bases into wildcards
			}
			my @matches = ();
			my %close = ();
			foreach my $barcode (keys(%CellBarcodes)) {
				if ($barcode =~/$CellID/) { # Match but with uncertainty
					push(@matches, $barcode);
				} else {
					if (scalar(@matches == 0)) { # Count mismatches
						my $count = ( $barcode ^ $CellID ) =~ tr/\0//;
						if ($count >= length($barcode)-$MAXmismatch) { # Allow upto 2 mismatches
							$close{$barcode} = $count;
						}
					}
				}
			}
			# If exact matches with uncertainty then give those priority, otherwise keep the most similar expected barcodes
			if (scalar(@matches) == 0 && scalar(keys(%close)) > 0) { # Has 1 or 2 mismatches
				my $max = my_max(values(%close)); # Closest match
				$mismatches = length($CellID)-$max;
				foreach my $code (keys(%close)) {
					if ($close{$code} == $max) {
						push(@matches,$code);
					}
				}
			}
			if (scalar(@matches) == 1) { # single best match
				$CellID = $matches[0];
				$Mismatch++;
			} elsif (scalar(@matches) > 1) { #More than one equally good match
				$AmbiguousCell++;
				next;
			} else { # No match
				$NotPossibleCell++;
				next;
			}
		} else { # Exact match
			$ExactMatch++;
		} 

		} #If known barcodes

		if ($ARGV[3] ne "UNKNOWN") {
			# Has Acceptable Barcode
			<$ifh1>;<$ifh2>; #+'s
			my $file1qual = <$ifh1>;
			my $file2qual = <$ifh2>;
			my $handle = $ofhs{$CellBarcodes{$CellID}};
			print $handle "$readname:$UMI\n$read+\n$file2qual";
			$OutputReads++;
		} else {
			$CellBarcodes{$CellID}++;
		}
	} else {next;}
}
if ($ARGV[3] ne "UNKNOWN") {
	print STDERR "
	Doesn't match any cell: $NotPossibleCell
	Ambiguous: $AmbiguousCell
	Exact Matches: $ExactMatch
	Contain mismatches: $Mismatch
	Input Reads: $total_reads
	Output Reads: $OutputReads\n";
	close($ifh1);
	close($ifh2);
	foreach my $ofh (keys(%ofhs)) {close($ofhs{$ofh});}
} else {
	print STDERR "Bad UMI: $BadUMI\n";
	my @Codes = sort { $CellBarcodes{$b} <=> $CellBarcodes{$a} } keys(%CellBarcodes);
	foreach my $code (@Codes) {
		print "$code ".$CellBarcodes{$code}."\n";
	}
}

sub my_max {
	if (scalar(@_) == 1) {return($_[0])};
	my $max = shift;
	foreach my $ele (@_) {
		if ($ele > $max) {$max = $ele;}
	}
	return($max);
}
