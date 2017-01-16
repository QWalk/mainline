
#get the actual polarization from two sets of files.  Uses
#get_zpol.pl; they have to be modified a bit, but hopefully 
#these will be improved after some time..

$au2Cm=57.216; 
$pi=3.14159;
$a=4.00/.529;
$c=4.036/.529;
$base1="centro";
$base2="dist";


$ion1=`grep "ionic polarization" ${base1}5.vmc.o | awk '{print \$5}'`;
$ion2=`grep "ionic polarization" ${base2}5.vmc.o | awk '{print \$5}'`;

chomp($ion1); chomp($ion2);
print "ionic contribution : $ion1  $ion2\n";

$iondiff=$ion2-$ion1;

@phase=split(' ', `perl get_zpol.pl | grep angle`);
$phasediff=$phase[1]; #$phase2[1]-$phase1[1];
$phaseerr=$phase[2]; #sqrt($phase2[3]**2+$phase1[3]**2);


$ne=1;

$vol=$a*$a*$c;
$area=$a*$a;


print "volume: $vol  area: $area\n";
$fac=$au2Cm*$ne/(2*$pi*$area);

$elec_contrib=$phasediff*$fac;
$ion_contrib=$au2Cm*$iondiff/($vol);
$err=$phaseerr*$fac;

$quantum=2*$pi*$fac;
print "quantum: $quantum \n";
print "electronic: $elec_contrib C/m^2  ionic: $ion_contrib  C/m^2\n";
$total=$ion_contrib-$elec_contrib;
print "difference: $total +/- $err C/m^2 \n";
