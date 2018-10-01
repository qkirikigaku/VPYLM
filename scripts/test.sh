#!bin/bash
export KMP_AFFINITY=disabled

for i in `seq 1 10`; do
    mkdir result/test_${i}
    ./VPYLM test Upstream $i
    ./VPYLM test Downstream $i
done
