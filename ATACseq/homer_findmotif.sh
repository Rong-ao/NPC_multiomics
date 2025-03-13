#!/bin/bash

lists=(H9_ESC NPC_D1 NPC_D2 NPC_D3 NPC_D6)
for list in ${lists[*]}
do
    echo ${list} begin
    findMotifsGenome.pl ${list}_merge/${list}_summits.bed hg38 ${list}_homer/ -size 200 -mask
    echo ${list} done
done