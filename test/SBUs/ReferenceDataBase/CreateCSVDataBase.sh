#!/bin/bash
cat *_si.dat | sed 's/,/ /g' | sed 's/  */ /g' | sed 's/^[[:space:]]*//g' | sed 's/[[:blank:]]*$//' | sed 's/ /,/g' > cosa.csv
