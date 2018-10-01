#!bin/bash

for i in `seq 1 99`; do
    python runvpylm.py ${i}
done
