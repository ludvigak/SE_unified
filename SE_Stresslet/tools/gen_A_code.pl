#!/usr/bin/perl

$s = '';
foreach $m (1..3) {
	foreach $j (1..3) {
		$s .= "A($j,$m)=";
		$s .= "C*x($j)*x($m)*xdotn";
		$s .= " + D*( x($m)*nvec($j) + x($j)*nvec($m) ";
		if ($m==$j) {
			$s .= "+ xdotn";
		}
		$s .= ");\n";
	}
}
$s =~ s/\++-/-/g;
print $s;
