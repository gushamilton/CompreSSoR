#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J sharebench -c 8 --mem=64G -t 01:00:00
#SBATCH -o /user/work/fh6520/showcase/sharebench/logs/bench-%x-%j.out
# GRID=small|big|confirm. LIB = the library under test. Each cell is a fresh Rscript; modes one/share interleaved in random order per replicate.
module load languages/R/4.5.1
grep -m1 'model name' /proc/cpuinfo | grep -q "Gold 6226R" || { echo "wrong CPU on $(hostname)"; exit 3; }
export OMP_NUM_THREADS=1
L=${LIB:?set LIB to the library under test}; B=/user/work/fh6520/showcase/sharebench/bench_share.R
MODES=${MODES:-"one share auto"}
case ${GRID:-small} in
  small) CELLS="showcase:20:1000:1 showcase:20:1000:8 showcase:20:100000:1 showcase:20:100000:8 showcase:50:1000:1 showcase:50:1000:8 showcase:50:100000:1 showcase:50:100000:8 mix:25:1000:8 mix:30:1000:8 mix:40:1000:8 ukbppp:20:6000:1 ukbppp:20:6000:8 ukbppp:60:6000:8 ukbppp:160:6000:8 ukbppp:300:6000:8";;
  big) CELLS="ukbppp:2940:6000:8";;
  confirm) CELLS="showcase:20:1000:8 showcase:20:100000:1 showcase:50:1000:1 mix:25:1000:8 ukbppp:20:6000:8 ukbppp:300:6000:8";;
esac
for rep in 1 2 3; do for c in $CELLS; do IFS=: read set n k t <<< "$c"
  for m in $(printf '%s\n' $MODES | shuf --random-source=<(yes $rep$n$k$t)); do
    /usr/bin/time -f "TIME,%e,%M" Rscript $B $L $set $n $k $t $m 2>&1 | grep -E "^RESULT|^TIME|Error" | paste -sd, | sed "s/^/rep$rep,/"
  done; done; done
