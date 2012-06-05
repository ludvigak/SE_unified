#!/usr/bin/perl

$s = '';
$s2 = '';
foreach $j (0..2) {
	foreach $m (0..2) {
		if ($j<=$m) {
		$s .= "A[$j][$m]=";
		$s .= "C*x[$j]*x[$m]*xdotn";
		$s .= " + D*( x[$m]*n[$j] + x[$j]*n[$m] ";
		if ($m==$j) {
			$s .= "+ xdotn";
		}
		$s .= ");\n";
		}
		else {
		$s2 .= "A[$j][$m]=A[$m][$j];\n";
		}
	}
}
$s .= $s2;

$s.= "-----------------------------\n";
foreach $k1 (0..2) {
	foreach $k2 (0..2) {
		$j = $k1;
		$m = $k2;
		if ($m<=$j) {
			$j = $k2;
			$m = $k1;
		}			
		$s .= "tmp[m+nidx*".($k1+3*$k2)."] = tmp[m+nidx*".($k1+3*$k2)."] + A[$j][$m];\n";
	}
}

$s =~ s/\++-/-/g;
print $s;

