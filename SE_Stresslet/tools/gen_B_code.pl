#!/usr/bin/perl

$s = '';
foreach $m (1..3) {
	foreach $l (1..3) {
		foreach $j (1..3) {
			$s .= "Bi($j,$l,$m)=";
			if($l>$j) {
				$s .= "Bi($l,$j,$m);\n";
			}
			elsif($m>$l) {
				$s .= "Bi($j,$m,$l);\n";
			}
			else {
				$s .= "psi*(";
				$s .= "k($m)+" if $j==$l;
				$s .= "k($j)+" if $m==$l;
				$s .= "k($l)+" if $j==$m;
				$s .= "-2*k($j)*k($l)*k($m)*k2i";
				$s .= ");\n";
			}
		}
	}
}
$s =~ s/\++-/-/g;
print $s;
