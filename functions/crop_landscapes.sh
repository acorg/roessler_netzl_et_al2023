#!/bin/bash

CURR_DIR=`pwd`

cd $CURR_DIR
cd ./figures/landscapes/gmt_landscapes

PNG_LIST=`ls -lart *.png | awk '{print $9}'`

echo "Cropping landscapes"

for f in $PNG_LIST; do
 #   convert ${f} -crop 2200x1200+500+600 ${f}
    convert ${f} -gravity center -crop 1950x1100-140+150 ${f}
done
