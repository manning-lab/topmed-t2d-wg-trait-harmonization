#!/bin/bash
# $1 = f.dir
# $2 = out.pref
# $3 = trait 
# $4 = continuous covariates for summary plot

# get a unique identified for the users current branch, this links back to the last commit on the branch that the user is working. Can be used to see the exact code state of the code (to the last remote commit)
git_tag=$(git describe --always --dirty)
user=$(whoami)
curdate=$(date +"%m.%d.%Y")

if [[ $git_tag == *dirty ]] ; then
	echo "You have uncommitted changes. Script will not run until changes are committed. Current Tag: $git_tag"
	exit 1
fi

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
echo "output file prefix : $2.$user.$curdate"
echo "Desired trait to removed duplicates : $3"
echo "\n"

echo "harmonization:\n" > $1/$2.$user.$curdate.harm.stdout.txt
# make the pooled trait file
(R --vanilla --args $1 $2.$user.$curdate < gly_traits_harmonization_021318.R >> $1/$2.$user.$curdate.harm.stdout.txt 2>&1) &
echo "Running harmonization" 
spinner $! 
echo "Done!"
echo "\n"

echo "\nduplicates:\n" >> $1/$2.$user.$curdate.harm.stdout.txt

(R --vanilla --args $1 $2.$user.$curdate TOPMEDID $2.$user.$curdate.csv $3 < duplicates_v2_TM_022218.R >> $1/$2.$user.$curdate.harm.stdout.txt 2>&1) &
echo "Running duplicates script" 
spinner $! 
echo "Done!"


(R --vanilla --args $1 $2.$user.$curdate TOPMEDID $2.$user.$curdate.no.duplicates.csv $3 clusters.list.sqrt_v4.RData < Freeze5b.FI_FG.PostProcessing.GitHub.R >> $1/$2.$user.$curdate.harm.stdout.txt 2>&1) &
echo "Running Post-Processing script" 
spinner $! 
echo "Done!"

(R --vanilla --args $2.$user.$curdate.for_analysis.csv $3 $4 sex $2.$user.$curdate.for_analysis topmed_project ancestry $1 < PhenotypeSummary.R T2D_FG>> $1/$2.$user.$curdate.harm.stdout.txt 2>&1) &
echo "Running Trait Summary script" 
spinner $! 
echo "Done!"

echo "git commit tag of current code : $git_tag" > $1/$2.$user.$curdate.log
echo "harmonized phenotypes with duplicates : $1/$2.$user.$curdate.csv" >> $1/$2.$user.$curdate.log
echo "harmonized phenotypes without duplicates : $1/$2.$user.$curdate.no.duplicates.csv" >> $1/$2.$user.$curdate.log
echo "list of NWD ids that were removed : $1/$2.$user.$curdate.removed.IDs.txt" >> $1/$2.$user.$curdate.log
echo "duplicates info table : $1/$2.$user.$curdate.duplicates.txt" >> $1/$2.$user.$curdate.log
echo "stdout and stderr file : $1/$2.$user.$curdate.harm.stdout.txt" >> $1/$2.$user.$curdate.log
