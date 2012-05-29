#!/usr/bin/perl

$s = '';
foreach $m (0..2) {
	foreach $j (0..2) {
		$s .= "A[$j][$m]=";
		$s .= "C*x[$j]*x[$m]*xdotn";
		$s .= " + D*( x[$m]*n[$j] + x[$j]*n[$m] ";
		if ($m==$j) {
			$s .= "+ xdotn";
		}
		$s .= ");\n";
	}
}
$s =~ s/\++-/-/g;
print $s;
