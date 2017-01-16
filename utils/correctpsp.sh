#!/bin/sh
if [[ $# > 0 ]]
then
    echo "Correcting $1"
    cp $1 $1.backup
    a=`grep "BASIS {" $1 | sed "s/BASIS {//g" | sed "N;s/\n//"`
    getpsp.py INPUT "$a" > psp
    sed -n  -n "1,/PSEUDO/p" $1.backup > $1 
    sed -i -e "/PSEUDO {/r psp" -e "//d" $1
else
    for f in `ls qwalk_*.sys`
    do
	echo "Correcting $f"
	cp $f $f.backup
	a=`grep "BASIS {" $f | sed "s/BASIS {//g" | sed "N;s/\n//"`
	getpsp.py INPUT "$a" > psp
	sed -n  -n "1,/PSEUDO/p" $f.backup > $f 
	sed -i -e "/PSEUDO {/r psp" -e "//d" $f
    done
fi
