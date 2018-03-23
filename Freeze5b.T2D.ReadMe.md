# T2D Harmonization readme

## Prereqs

1. Put all of the study-specific phenotype files, below, in a single directory. Note the absolutely path to that directory.

2. Make sure your github repo is up-to-date and that you're on the correct branch.

```bash
cd /my/local/git/repo/path/topmed-traitHarmonization/ 
git stash 
git fetch
git checkout master
git pull origin master
```

3. Make sure that youre in the local git repo directory (something like /my/compy/topmed-traitHarmonization/)

## How to Run

### From command line

inputs:

1. /my/pheno/dir/ : absolute path to study specific phenotype files (from step 1 in prereqs)
2. prefix : desired prefix for output file
3. trait : if you'd like to removed duplicates based on trait values, input the column name in the phenotype file for the trait of interest


run:
```bash
sh run_t2d_harmonization.sh /my/pheno/dir/from/step/1/ prefix trait

```

ouputs:
1. /my/pheno/dir/prefix.yourusername.MM.DD.YYYY.csv : harmonized phenotype file with duplicates
2. /my/pheno/dir/prefix.yourusername.MM.DD.YYYY.duplicates.txt : table of phenotype data for duplicates
3. /my/pheno/dir/prefix.yourusername.MM.DD.YYYY.removed.IDs.txt : IDS and reason for removal of duplicates
4. /my/pheno/dir/prefix.yourusername.MM.DD.YYYY.no.duplicates.csv : harmonized phenotype file without duplicates
5. /my/pheno/dir/harm.stdout : stdout and stderr of both harmonization and duplicate removal
6. /my/pheno/dir/prefix.yourusername.MM.DD.YYYY.log : logfile with filepaths to outputs and unique identified for code used