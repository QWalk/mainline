
#use strict;
sub get_realimag;
sub magnitude;

#We calculate the angle difference between two Z=<exp(i*G*X)>
#as cos(theta)=a*b/ab.  This seems to be the best way to do it,
#since calculating the angles directly seems to be very noisy..
#in each array, [0] is the average value, [1] is the error bar, and
#[2] is the bias, if any
$string="z_pol2";
$prefix=$ARGV[0];
$suffix=".vmc.log";

($areal[0],$areal[1],$aimag[0], $aimag[1])=get_realimag("centro", $suffix);
print "centro real/imag $areal[0]  $aimag[0] errors  $areal[1]  $aimag[1] \n";
@a=magnitude(@areal, @aimag);
print "centro magnitude $a[0]  $a[1]  bias $a[2] \n\n";

($breal[0],$breal[1],$bimag[0], $bimag[1])=get_realimag("dist", $suffix);
print "dist real/imag $breal[0]  $bimag[0] errors  $breal[1]  $bimag[1] \n";
@b=magnitude(@breal, @bimag);
print "dist magnitude $b[0]  $b[1]  bias $b[2] \n\n";

$ab[0]=$a[0]*$b[0];
$ab[1]=abs($a[0])*$b[1]+abs($b[0])*$a[1];
$ab[2]=0;

$adotb[0]=$areal[0]*$breal[0]+$aimag[0]*$bimag[0];
$adotb[1]=$areal[1]*$breal[0]+$areal[0]*$breal[1]
    +$aimag[1]*$bimag[0]+$aimag[0]*$bimag[1];
$adotb[2]=0;

@cos_2=divide(@adotb,@ab);
print "\n\n  adotb @adotb   ab  @ab \n";
print "cos method 2 @cos_2 \n";

@angle=acos(@cos_2);
print "\n angle @angle\n";


sub acos {
   my($x) = shift(@_);
   my($xerr)=shift(@_);
   my $ret = atan2(sqrt(1 - $x**2), $x);
   my $err=$xerr*1/(sqrt(1-$x**2));
   my $bias=0;
   return ($ret,$err,$bias);
}

sub divide {
    my @x,@y,@z;
    for(my $i=0; $i < 3; $i++) {
	$x[$i]=shift(@_);
    }
    for(my $i=0; $i < 3; $i++) {
	$y[$i]=shift(@_);
    }
    $z[0]=$x[0]/$y[0];
    $z[1]=$x[1]/abs($y[0])+$y[1]*abs($z[0]);
    $z[2]=$y[1]*$y[1]*2.0*$x[0]/$y[0]**3;
    return @z;
}

sub magnitude { 
    my @real,@imag,$mag, $magerr, $magbias;
    $real[0]=shift(@_);
    $real[1]=shift(@_);
    $imag[0]=shift(@_);
    $imag[1]=shift(@_);
    $mag=sqrt($real[0]*$real[0]+$imag[0]*$imag[0]);
    $magerr=(abs($real[0])*$real[1]+abs($imag[0])*$imag[1])/$mag;
    $magbias=.5*($real[1]**2*($mag**2-$real[0]**2)
		 +$imag[1]**2*($mag**2-$imag[0]**2))/$mag**3;
    return ($mag, $magerr, $magbias);
}

sub get_realimag {
    my @gamma, @k100, @k110, @k001, @k101, @lpoint;
    my @real, @imag;
    my $prefix=shift(@_);
    my $suffix=shift(@_);
    @gamma=split(' ',`gosling ${prefix}0${suffix} | grep z_pol2`);
    @k100=split(' ',`gosling ${prefix}1${suffix} | grep z_pol2`);
    @k110=split(' ',`gosling ${prefix}2${suffix} | grep z_pol2`);
    @k001=split(' ',`gosling ${prefix}3${suffix} | grep z_pol2`);
    @k101=split(' ',`gosling ${prefix}4${suffix} | grep z_pol2`);
    @lpoint=split(' ',`gosling ${prefix}5${suffix} | grep z_pol2`);
    
    $real[0]=(.125*$gamma[1]
	      +.250*$k100[1]
	      +.125*$k110[1]
	      +.125*$k001[1]
	      +.25*$k101[1]
	      +.125*$lpoint[1]);
    $real[1]=sqrt(.125*.125*$gamma[6]**2
		  +.250*.250*$k100[6]**2
		  +.125*.125*$k110[6]**2
		  +.125*.125*$k001[6]**2
		  +.25*.25*$k101[6]**2
		  +.125*.125*$lpoint[6]**2);
    
    $imag[0]=(.125*$gamma[4]
	      +.250*$k100[4]
	      +.125*$k110[4]
	      +.125*$k001[4]
	      +.25*$k101[4]
	      +.125*$lpoint[4]);
    
    $imag[1]=sqrt(.125*.125*$gamma[9]**2
		  +.250*.250*$k100[9]**2
		  +.125*.125*$k110[9]**2
		  +.125*.125*$k001[9]**2
		  +.25*.25*$k101[9]**2
		  +.125*.125*$lpoint[9]**2);
    return (@real, @imag);
}
