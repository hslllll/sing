#!/usr/bin/env bash
set -e

echo "cd `pwd`; bash scripts/benchmark.sh y" | qsub -N benchmark-y -l nodes=node01:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-y.log
echo "cd `pwd`; bash scripts/benchmark.sh a" | qsub -N benchmark-a -l nodes=node01:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-a.log
echo "cd `pwd`; bash scripts/benchmark.sh m" | qsub -N benchmark-m -l nodes=node01:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-m.log
echo "cd `pwd`; bash scripts/benchmark.sh h" | qsub -N benchmark-h -l nodes=node01:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-h.log
echo "cd `pwd`; bash scripts/benchmark-gatk.sh y" | qsub -N benchmark-gy -l nodes=node03:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-gy.log
echo "cd `pwd`; bash scripts/benchmark-gatk.sh a" | qsub -N benchmark-ga -l nodes=node03:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-ga.log
echo "cd `pwd`; bash scripts/benchmark-gatk.sh m" | qsub -N benchmark-gm -l nodes=node03:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-gm.log
echo "cd `pwd`; bash scripts/benchmark-gatk.sh h" | qsub -N benchmark-gh -l nodes=node03:ppn=16 -v WORKDIR=`pwd` -j oe -o ~/log/benchmark-gh.log

echo "ALL JOBS SUMBMITTED!"
