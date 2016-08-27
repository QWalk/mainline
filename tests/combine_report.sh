#!/bin/bash
cat h/report.csv > report.csv
tail -n +2 h2/report.csv >> report.csv
tail -n +2 n2/report.csv >> report.csv
