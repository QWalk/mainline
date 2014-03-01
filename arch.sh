#!/bin/sh
WORKDIR=$HOME/Documents/Research/vo2/Programs/qwalk-zhh
function compileconverter() {
    cd $WORKDIR/src/converter
    make 
    cd -
}

function runqwalk()
{
    if [[ $# < 2 ]] 
    then
	$WORKDIR/src/qwalk-Linux-mpi $1
	echo "$WORKDIR/src/qwalk-Linux-mpi $1"
    else
	echo "mpirun -np $2 $WORKDIR/src/qwalk-Linux-mpi $1"
	mpirun -np $2 $WORKDIR/src/qwalk-Linux-mpi $1
    fi
}

function runqwalk_old()
{
    if [[ $# < 2 ]] 
    then
	qwalk $1
    else
	mpirun -np $2 qwalk $1
    fi
}


function rungosling()
{
    $WORKDIR/src/gosling-Linux-mpi $@
}

function runcrystal2qmc() {
    $WORKDIR/src/converter/crystal2qmc-Linux-mpi $@
}

function runcrystal_post() {
    sed -f ~/Templates/crystalout.sed $@
}