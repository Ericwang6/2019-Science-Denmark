#!/bin/bash
dir=(24 48 72 96 144 192 240 288 336)
for((i=0;i<${#dir[*]};i++))
do
    d=${dir[i]};
    ./ModelSelection.py $d
    echo $d" Finished."
done
