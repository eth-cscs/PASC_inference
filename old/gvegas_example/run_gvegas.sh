#!/bin/sh

ncall=e06n1
itmx=5
blk=64

if [ $# -ge 1 ];then
    ncall=$1
fi
if [ $# -ge 2 ];then
    itmx=$2
fi
if [ $# -ge 3 ];then
    blk=$3
fi

logfile=log_gvegas_example

date=`date +%Y%m%d%H%M%S`
if [ -f ${logfile} ]; then
    mv ${logfile} ${logfile}_${date}
fi

gvegas_example -ncall=${ncall} -itmx=${itmx} -blk=${blk} 1> ${logfile} 2>&1

