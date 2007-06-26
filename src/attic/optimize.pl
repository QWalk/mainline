#!/bin/perl

$basename="test";

$nconfig=10;
$eref=-1052.6;
$nregen=3;

$niterations=300;
$scrdir="/var/tmp";

$nblock=5;
$nstep=10;
$ndecorr=2;
$timestep=1.0;

$configdir="config";
$configname="$configdir/$basename.config";
system("mkdir $configdir");


$sysname=$basename . ".sys";
$jastname=$basename . ".jast";
$slatname=$basename . ".slater";

open(hfout, ">$basename.hf");
print hfout
"METHOD {
  VMC
  NBLOCK   $nblock
  NSTEP    $nstep
  NDECORR   $ndecorr
  TIMESTEP $timestep
  NCONFIG  $nconfig
  STORECONFIG  $configname
  EREF $eref
}



INCLUDE  $sysname
TRIALFUNC {
 INCLUDE $slatname
}
";
close(hfout);  

for($i=0; $i < $nregen; $i++) {
    
    open(optout, ">$basename.opt$i");
    print optout 
"METHOD {
  OPTIMIZE
  ITERATIONS $niterations
  READCONFIG $configname
  NCONFIG  $nconfig
  EREF $eref
  MINFUNCTION VARIANCE
  PSEUDOTEMP $scrdir/$basename.pseudo
}

INCLUDE $sysname  
";
    if($i==0) {
	print optout "TRIALFUNC { \n SLATER-JASTROW \n";
	print optout "  WF1 { INCLUDE $slatname } \n";
	print optout "  WF2 { INCLUDE $jastname } \n";
	print optout "}\n";
    }
    else {
	$i2=$i-1;
	print optout "TRIALFUNC { INCLUDE $basename.opt$i2.wfout }\n";
    }

    close(optout);

    open(vmcout, ">$basename.vmc$i");
    print vmcout 
"
METHOD {
  VMC
  NBLOCK   $nblock
  NSTEP    $nstep
  NDECORR   $ndecorr
  TIMESTEP $timestep
  NCONFIG  $nconfig
  STORECONFIG  $configname

  #uncomment the following to read a configuration
  READCONFIG  $configname
  EREF $eref
}



INCLUDE $sysname
TRIALFUNC {
INCLUDE $basename.opt$i.wfout
}
";
    close(vmcout);
}

