#!/bin/bash
# $1 = f.dir
# $2 = source.file
# $3 = out.pref
# $4 = trait 

spinner()
{
    local pid=$1
    local delay=0.75
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

echo "phenotype file directory : $1"
echo "filepath script path : $2"
echo "output file prefix : $3"
echo "Desired trait to removed duplicates : $4"
echo "Output harmonized file : $3.csv"
echo "\n"

# make the pooled trait file
(R --vanilla --args $1 $2 $3 < Harmonization.19JAN2017.GitHub.R > /dev/null 2>&1) &
echo "Running harmonization" 
spinner $! 
echo "Done!"
echo "\n"

(R --vanilla --args $1 topmedid $3.csv $4 $3 < duplicates_v2_TM_022218.R > /dev/null 2>&1) &
echo "Running duplicates script" 
spinner $! 
echo "Done!"