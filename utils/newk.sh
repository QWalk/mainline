#!/bin/sh 
if [ $# -lt 3 ]; then
    echo "Syntax: newk.sh crystal_out newk_out output"
    exit 1
fi
sed -n -e '1,/FINAL EIGENVECTORS/p' -e 's/FINAL EIGENVECTORS/HAMILTONIAN EIGENVECTORS/g' $1 > $3
sed -n -e '/HAMILTONIAN EIGENVECTORS/, $p' $2 >> $3
