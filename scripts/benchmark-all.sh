#!/usr/bin/env bash
set -e

echo "cd \`pwd\`; python3 scripts/benchmark.py 8 all y" | qsub -N benchmark-y -l nodes=node04:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-y.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 all a" | qsub -N benchmark-a -l nodes=node04:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-a.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 all m" | qsub -N benchmark-m -l nodes=node04:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-m.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 all b" | qsub -N benchmark-b -l nodes=node04:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-b.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 all h" | qsub -N benchmark-h -l nodes=node04:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-h.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 gatk y" | qsub -N benchmark-gy -l nodes=node03:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-gy.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 gatk a" | qsub -N benchmark-ga -l nodes=node03:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-ga.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 gatk m" | qsub -N benchmark-gm -l nodes=node03:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-gm.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 gatk b" | qsub -N benchmark-gb -l nodes=node03:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-gb.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 gatk h" | qsub -N benchmark-gh -l nodes=node03:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-gh.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 minimap y" | qsub -N benchmark-my -l nodes=node02:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-my.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 minimap a" | qsub -N benchmark-ma -l nodes=node02:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-ma.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 minimap m" | qsub -N benchmark-mm -l nodes=node02:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-mm.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 minimap b" | qsub -N benchmark-mb -l nodes=node02:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-mb.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 minimap h" | qsub -N benchmark-mh -l nodes=node02:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-mh.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 strobe y" | qsub -N benchmark-sy -l nodes=node01:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-sy.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 strobe a" | qsub -N benchmark-sa -l nodes=node01:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-sa.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 strobe m" | qsub -N benchmark-sm -l nodes=node01:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-sm.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 strobe b" | qsub -N benchmark-sb -l nodes=node01:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-sb.log
echo "cd \`pwd\`; python3 scripts/benchmark.py 8 strobe h" | qsub -N benchmark-sh -l nodes=node01:ppn=8 -v WORKDIR=\`pwd\` -j oe -o ~/log/benchmark-sh.log

echo "ALL JOBS SUMBMITTED!"
