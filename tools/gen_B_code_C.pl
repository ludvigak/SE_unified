#!/usr/bin/perl

$s = '';
foreach $j (0..2) {
	foreach $l (0..2) {
		foreach $m (0..2) {
			$s .= "Bi[".($j*9+$l*3+$m)."]="; # Vector indexing
#			if($l>$j) {
#				$s .= "Bi($l,$j,$m);\n";
#			}
#			elsif($m>$l) {
#				$s .= "Bi($j,$m,$l);\n";
#			}
#			else {
				$s .= "psi*(";
				$s .= "k[$m]+" if $j==$l;
				$s .= "k[$j]+" if $m==$l;
				$s .= "k[$l]+" if $j==$m;
				$s .= "-2*k[$j]*k[$l]*k[$m]*k2i";
				$s .= ");\n";
#			}
		}
	}
}
$s =~ s/\++-/-/g;
print $s;
