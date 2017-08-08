#!/bin/sh

nreg=64
if [ $# -ge 1 ]; then
    nreg=$1
fi

make clean;
make NREG=${nreg} >& makelog
