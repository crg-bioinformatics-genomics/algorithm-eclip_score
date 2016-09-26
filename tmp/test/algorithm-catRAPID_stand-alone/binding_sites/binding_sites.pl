use warnings;
use strict;
use Data::Dumper;
use Statistics::Basic qw(:all);

my $file;

### RBDs database

my @RBDs;
$file = 'RBDs/Gerstberger2014.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($name)=$_;
	push @RBDs, $name unless grep{$name eq $_}@RBDs;
}
close FILE;

$file = 'RBDs/NPIDB.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($name,$type)=(split/\s+/,$_)[0,2];
	next if $type eq "dna";
	push @RBDs, $name unless grep{$name eq $_}@RBDs;
}
close FILE;

my %CirilloSeq;
$file = 'RBDs/Cirillo2016.seq';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($prot,$seq)=split/\s+/,$_;
	$CirilloSeq{$prot} = $seq;
}
close FILE;

my %CirilloCoord;
$file = 'RBDs/Cirillo2016.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($prot,$coord)=(split/\s+/,$_)[0,2];
	push @{$CirilloCoord{$prot}}, $coord;
}
close FILE;

### RBDs detection in input protein

my %prot;
$file = '../../protein/outfile';
my @additionalRBDs;
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($id,$seq)=split/\s+/,$_;
	$prot{$id} = $seq;
	foreach my $p(keys %CirilloSeq){
		if($seq eq $CirilloSeq{$p}){
			foreach my $d(@{$CirilloCoord{$p}}){
				push @additionalRBDs, $d;
			}
		}else{
			next;
		}
	}
}
close FILE;

open OUT, '>', '../../tmp/prot.tmp' or die $!;
foreach my $id (keys %prot){
	print OUT join("\n",">".$id,$prot{$id}),"\n";
}
close OUT;

`hmmscan --cut_ga --noali --domtblout ../../tmp/pfam.out ../../../../../Pfam/Pfam-A.hmm ../../tmp/prot.tmp &> /dev/null`;

$file = '../../tmp/pfam.out';
my %pfam;
open PFAM, $file or die $!;
while(<PFAM>){
	next if $. < 4;
	if($_!~m/^\#/g){
		my @fields = split(" ",$_);
		if(grep{$fields[0] eq $_}@RBDs){
			my $coord = join("-",$fields[17],$fields[18]);
			$pfam{$coord} = $fields[0];
		}else{
			next;
		}
	}
}

foreach my $coord(@additionalRBDs){
	$pfam{$coord} = "Cirillo2016" unless grep{$coord eq $_}keys %pfam;
}

### selection of fragments based on RBD coverage

my @fragms;
$file = '../protein_library_generator/prot.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($id)=(split/\s+/,$_)[0];
	my $coord;
	if($id eq "1-50"){
		$coord = "1-50";
	}else{
		($coord) = $id=~ m/_(\d+\-\d+)$/g;
	}
	push @fragms, $coord unless grep{$coord eq $_}@fragms;
}
close FILE;

my @selectedFragms;
foreach my $fragmCoord (@fragms){
	my($sp,$ep) = split/\-/,$fragmCoord;
	foreach my $rbdCoord (keys %pfam){
		my($sd,$ed) = split/\-/,$rbdCoord;
			my $out = relativeCov($sp,$ep,$sd,$ed);
			if((split/\t/,$out)[2] > 0.25){
				push @selectedFragms, (split/\t/,$out)[0] unless grep{(split/\t/,$out)[0] eq $_}@selectedFragms;
			}
	}
}


### create table for binding sites detection

my %pred;
$file = '../combine_parallel/pre-compiled/out.merged.posi.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($prot,$rna,$score)=(split/\s+/,$_)[0..2];
	my $protName;
	my $protCoord;
	my $rnaName;
	my $rnaCoord;
	if($prot ne "1-50"){
		($protName,$protCoord) = $prot =~ m/^(\S+)\_(\d+\-\d+)$/g;
	
	}else{
		$protName = "prot";
		$protCoord = $prot;
	}
	if($rna ne "1-50"){
		($rnaName,$rnaCoord) = $rna =~ m/^(\S+)\_(\d+\-\d+)$/g;
	}else{
		$rnaName = "rna";
		$rnaCoord = $rna;	
	}
	$pred{$protCoord}{$rnaCoord} = $score;
}
close FILE;

my %negatmp;
$file = '../combine_parallel/pre-compiled/out.merged.nega.txt';
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($prot,$rna,$score)=(split/\s+/,$_)[0..2];
	my $rnaName;
	my $rnaCoord;
	if($rna ne "1-50"){
                ($rnaName,$rnaCoord) = $rna =~ m/^(\S+)\_(\d+\-\d+)$/g;
        }else{
                $rnaName = "rna";
                $rnaCoord = $rna;
        }
	push @{$negatmp{$rnaCoord}}, $score;
}
close FILE;
my %nega;
foreach my $rnaCoord(keys %negatmp){
	$nega{$rnaCoord} = join("\t",mean(@{$negatmp{$rnaCoord}}),stddev(@{$negatmp{$rnaCoord}}));
}

my %table;
my @tmp;
foreach my $protCoord (keys %pred){
	foreach my $rnaCoord (keys %{$pred{$protCoord}}){
		my $rbd = 0;
		$rbd = 1 if grep{$protCoord eq $_}@selectedFragms;
		my $score = $pred{$protCoord}{$rnaCoord};
		my($mean,$stddev) = split/\s+/,$nega{$rnaCoord};
		my $norm = ($score-$mean)/$stddev;
		my $p = $protCoord;
		$p =~ s/\-/\t/g;
		my $r = $rnaCoord;
		$r =~ s/\-/\t/g;
		$table{$protCoord}{$rnaCoord} = join("\t",$p,$r,$score,$nega{$rnaCoord},$rbd,$norm);
		push @tmp, $norm if $rbd eq 1;
	}
}

my $ref_perc = percentile(\@tmp);
my %perc = %{$ref_perc};

print join("\t",qw/startProt endProt startRNA endRNA score negaMean negaStddv RBDs normScore percentileRBDs/),"\n";

foreach my $protCoord (keys %table){
	foreach my $rnaCoord (keys %{$table{$protCoord}}){
		my $norm = (split/\t/,$table{$protCoord}{$rnaCoord})[-1];
		if(defined $perc{$norm}){
			print join("\t",$table{$protCoord}{$rnaCoord},$perc{$norm}),"\n";
		}else{
			print join("\t",$table{$protCoord}{$rnaCoord},0),"\n";
		}
	}
}

###

sub relativeCov {

my $sp=shift; # fragm coord
my $ep=shift;
my $sd=shift; # RBD coord
my $ed=shift;

my $cov;
if($ep<$sd or $sp>$ed){
	$cov = 0;
}
if($ep<=$ed or $sp<=$sd){
	$cov = ($ep-$sd+1)/($ed-$sd+1);
}
if($ep>=$ed or $sp>=$sd){
	$cov = ($ed-$sp+1)/($ed-$sd+1);
}

if($cov>1){$cov = 1;}
if($cov<0){$cov = 0;}

return join("\t",$sp."-".$ep,$sd."-".$ed,$cov);

}

sub percentile {
	my %percentile;
	my %count;
	my @data = @{$_[0]};
	foreach my $datum (@data) {
    	++$count{$datum};
    }
	my $total = 0;
	foreach my $datum (sort { $a <=> $b } keys %count) {
    	$total += $count{$datum};
    	$percentile{$datum} = $total / @data;
	}
	return \%percentile;
}

