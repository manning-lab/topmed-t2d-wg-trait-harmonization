#!/bin/bash
# $1 = f.dir
# $2 = source.file
# $3 = out.pref
# $4 = trait 

# get a unique identified for the users current branch, this links back to the last commit on the branch that the user is working. Can be used to see the exact code state of the code (to the last remote commit)
git_tag=$(git rev-list HEAD | head -1)
user=$(whoami)
curdate=$(date +"%m.%d.%Y")

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
echo "\n"

echo "harmonization:\n" > $1/harm.stdout
# make the pooled trait file
(R --vanilla --args $1 $2 $3.$user.$curdate < Harmonization.19JAN2017.GitHub.R >> $1/harm.stdout 2>&1) &
echo "Running harmonization" 
spinner $! 
echo "Done!"
echo "\n"

echo "\nduplicates:\n" >> $1/harm.stdout

(R --vanilla --args $1 $3.$user.$curdate topmedid $3.$user.$curdate.csv $4 < duplicates_v2_TM_022218.R >> $1/harm.stdout 2>&1) &
echo "Running duplicates script" 
spinner $! 
echo "Done!"

echo "git commit tag of current code : $git_tag" > $1/$3.$user.$curdate.log
echo "harmonized phenotypes with duplicates : $1/$3.$user.$curdate.csv" >> $1/$3.$user.$curdate.log
echo "harmonized phenotypes without duplicates : $1/$3.$user.$curdate.no.duplicates.csv" >> $1/$3.$user.$curdate.log
echo "list of NWD ids that were removed : $1/$3.$user.$curdate.removed.IDs.txt" >> $1/$3.$user.$curdate.log
echo "duplicates info table : $1/$3.$user.$curdate.duplicates.txt" >> $1/$3.$user.$curdate.log
echo "stdout and stderr file : $1/harm.stdout" >> $1/$3.$user.$curdate.log