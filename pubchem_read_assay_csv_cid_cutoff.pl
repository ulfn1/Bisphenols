#!/usr/bin/perl


# MIT License

# Copyright (c) 2022 Ulf Norinder

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


my $unik = shift @ARGV;
my $unik2 = shift @ARGV;
my $unik3 = shift @ARGV;
my $unik4 = shift @ARGV;
my $unik5 = shift @ARGV;

if ($unik eq "") {
        print "You must give a pubchem comma sep input file name\n";
        exit;
}

if ($unik2 eq "") {
        print "You must give a pubchem tab sep cid smiles or sid smiles file name (smiles cid/sid)\n";
        exit;
}

if ($unik3 eq "") {
        print "You must indcate cutoff level for activity (0-1)\n";
        exit;
}

if ($unik4 eq "") {
        print "You must indcate smiles with sid (s) or cid (c)\n";
        exit;
}

$obj = 0;

open (IN, "<$unik2") || die "$0 cannot open Inputfile: $unik2\n";
while ($inline = <IN>) {
	chomp($inline);
	@arrayOneLine = split(/\t/,$inline);
	if($arrayOneLine[0] ne "") {
		if ($arrayOneLine[1] ne "") {
			$smi{$arrayOneLine[1],1} = $arrayOneLine[0];
		}
	}
}
close IN;

open (IN, "<$unik") || die "$0 cannot open Inputfile: $unik\n";
$inline = <IN>;
chomp($inline);
@arrayOneLine = split(/\,/,$inline);
$sizeoffile = @arrayOneLine;

$cid = -1;
$act = -1;
for ($abc = 0; $abc <$sizeoffile; $abc++) {
	if($unik4 eq "c") {
		if($arrayOneLine[$abc] eq "PUBCHEM_CID") {
			$cid = $abc;
		}
	}
	if($unik4 eq "s") {
		if($arrayOneLine[$abc] eq "PUBCHEM_SID") {
			$cid = $abc;
		}
	}
	if($arrayOneLine[$abc] eq "PUBCHEM_ACTIVITY_OUTCOME") {
		$act = $abc;
	}
}
if($cid < 0) {
	print "PUBCHEM_CID/SID field missing. Exiting\n";
        exit;
}
if($act < 0) {
	print "PUBCHEM_ACTIVITY_OUTCOME field missing. Exiting\n";
        exit;
}
open (OUT, ">$unik.data") || die "$0 cannot open Outputfile: $unik.data\n";
print OUT "Name\ttarget\tsmiles\n";
open (OUT2, ">$unik.smi") || die "$0 cannot open Outputfile: $unik.smi\n";
open (OUT3, ">$unik.target") || die "$0 cannot open Outputfile: $unik.target\n";
print OUT3 "Name\ttarget\n";

@cmpds = ();
@smiles= ();
while ($inline = <IN>) {
	chomp($inline);
	@arrayOneLine = split(/\,/,$inline);
	if($arrayOneLine[$cid] ne "") {
		push (@cmpds , $arrayOneLine[$cid]);
		if ($arrayOneLine[$act] ne "") {
			$aa = $smi{$arrayOneLine[$cid],1};
			$bb = $arrayOneLine[$cid];
			if ($aa eq "") {
				print "$unik, smiles for cid $arrayOneLine[$cid] is missing\n";
			} else {
				push (@smiles , $aa);
				if ($arrayOneLine[$act] eq "Active") {
					$data{$aa,1} = $data{$aa,1} + 1;
				}
				if ($arrayOneLine[$act] eq "Inconclusive") {
					$data{$aa,3} = $data{$aa,3} + 1;
				}
				if ($arrayOneLine[$act] eq "Inactive") {
					$data{$aa,2} = $data{$aa,2} + 1;
				}
				if ($bb ne "") {
					$smi{$aa,2} = $bb;	
				} else {
					print "$unik, cid $arrayOneLine[$cid] is missing\n";
				}
			}
		}
	}
}

my @unique = do { my %seen; grep { !$seen{$_}++ } @smiles };
close IN;

$activetot = 0;
$inactivetot = 0;
foreach (@unique) {
	$aa1 = $data{$_,1};
	if ($aa1 eq "") {
		$aa1 = 0;
	}
	$aa2 = $data{$_,2};
	if ($aa2 eq "") {
		$aa2 = 0;
	}
	$aa3 = $data{$_,3};
	if ($aa3 eq "") {
		$aa3 = 0;
	}
	print "$smi{$_,2}, $_, $aa1, $aa2, $aa3\n";
	$tot = $aa1 + $aa2 + $aa3;
	$active = $aa1/$tot;
	$inactive = $aa2/$tot;
	if ($aa1 != 0 and $active > $unik3) {
		$activetot = $activetot + 1;
		print OUT "$smi{$_,2}\t1\t$_\n";
		print OUT2 "$_\t$smi{$_,2}\n";
		print OUT3 "$smi{$_,2}\t1\n";
	}
	if ($aa2 != 0 and $inactive > $unik3) {
		$inactivetot = $inactivetot + 1;
		print OUT "$smi{$_,2}\t-1\t$_\n";
		print OUT2 "$_\t$smi{$_,2}\n";
		print OUT3 "$smi{$_,2}\t-1\n";
	}
}


close OUT;
close OUT2;
close OUT3;

open (OUT, ">$unik.cls") || die "$0 cannot open Outputfile: $unik.cls\n";
print "$unik, Actives: $activetot , Inactives: $inactivetot\n";
print OUT "$unik, Actives: $activetot , Inactives: $inactivetot\n";
close OUT;


exit;
