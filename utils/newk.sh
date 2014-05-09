#!/bin/sh 
if [ $# -lt 3 ]; then
    echo "Syntax: newk.sh crystal_out newk_out output"
    exit 1
fi
head -50000 $1 > $1.tmp
sed -n -e '1,/FINAL EIGENVECTORS/p' $1.tmp | sed 's/FINAL EIGENVECTORS/NEWK EIGENVECTORS/g' > $3
sed -n -e '/HAMILTONIAN EIGENVECTORS/, $p' $2 >> $3
