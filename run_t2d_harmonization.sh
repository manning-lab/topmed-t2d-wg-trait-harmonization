#!/bin/bash
# $1 = f.dir
# $2 = source.file
# $3 = out.pref

# make the pooled trait file
R --vanilla --args $1 $2 $3 < Harmonization.19JAN2017.GitHub.R

# removed the duplicates

Rscript duplicates.R