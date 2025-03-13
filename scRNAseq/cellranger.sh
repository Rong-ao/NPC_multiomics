#!/bin/bash

lists=(H9_ESC NPC_24h NPC_44h NPC_46h NPC_48h NPC_50h NPC_52h NPC_54h NPC_72h NPC_D6)
for list in ${lists[*]}
do
    echo ${list} begin
    cellranger-arc count --id=${list} \
                        --reference=/ref_data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                        --libraries=/data/${list}_library.csv \  # this is for multi-omics data
                        --localcores=16 \
                        --localmem=64
    echo ${list} done
done