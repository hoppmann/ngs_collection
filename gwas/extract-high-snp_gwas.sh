#!/bin/bash

head -n 1 $1 > $1.csv

awk '{ if ($4 >= 0.05 && $4 <= 0.95 && $15 > 1 ) print}' $1 | sort -t$'\t' -g -k10 | head -n 20 >> $1.csv
