#!/bin/bash

CURR_DIR=`pwd`

cd $CURR_DIR
cd ./figures/landscapes/gmt_landscapes

HTML_LIST=`ls -lart *.html | awk '{print $9}'`

echo "Capturing screenshots of html landscapes"

for f in $HTML_LIST; do
    new_name=${f%.*}.png
    echo $new_name
    open $f
    sleep 1
   # screencapture -x -w $new_name
    screencapture -x -m $new_name
done
