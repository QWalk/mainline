use File::Copy;
$installpath=$ARGV[0];
$platform=$ARGV[1];

@files=split(' ', `echo *-$platform converter/*-$platform`);
foreach $file (@files) { 
  $stripped=$file;
  $stripped=~ s/converter\///g;
  $stripped=~ s/-$platform//;
  system("mkdir -p $installpath");
  system("cp -fa $file  $installpath/$stripped");
}

