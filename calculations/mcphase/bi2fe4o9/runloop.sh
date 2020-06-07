#!/bin/bash 

mkdir -p results
cif2mcphas.pl bi2fe4o9.cif -pc $1 >& results/cif2mcphas.out
for ff in Fe*sipf
do
    sed "s/^L00/#L00/" $ff > t.s && echo >> t.s && echo truncate_matrix=0.025 >> t.s && mv t.s $ff
    ic1ion $ff && mv results/ic1ion.out results/ic1ion.out-$ff
    anisotropy_fibo 1 1  100000 -r $ff 0 0 0 0 0 0 >&/dev/null \
    && python get_aniso.py results/anisotropy.out >& results/anisodir-$ff \
    && mv results/anisotropy.out results/anisotropy.out-$ff
done >&/dev/null
for ff in results/ic1ion.out-Fe*; do echo $ff $(grep "^[0-9]" $ff | head -n6 | awk '{print $1}'); done > results/cef_en.out
for ff in ini tst j; do cp mcphas.${ff}-0 mcphas.${ff}; done
mcphasit -v >& results/mcphasit.out
mv results results-pc-$1-110
mkdir results
mcphasit -v >& results/mcphasit.out
mv results results-pc-$1-001
cp results-pc-$1-110/mcphas.fum ../bfo_110_mcphas.fum
cp results-pc-$1-001/mcphas.fum ../bfo_001_mcphas.fum
