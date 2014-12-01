#!/bin/bash

for((m=4;m<=12;m+=2));
do
python rewriteWMsmall.py 17 $m;
qsub runWMsmall.sh;
done

for((m=4;m<=12;m+=2));
do
python rewriteWMsmall.py 30 $m;
qsub runWMsmall.sh;
done

