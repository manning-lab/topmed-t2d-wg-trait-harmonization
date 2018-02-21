# T2D Harmonization readme

## Prereq's

1. Put all of the study-specific phenotype files, below, in a single directory. Note the absolutely path to that directory.

2. Make sure your github repo is up-to-date and that you're on the correct branch.

```bash
cd /my/local/git/repo/path/topmed-t2d-glycemia-public/ 
git stash 
git fetch
git checkout glycemic_traits_harmonization
git pull origin glycemic_traits_harmonization
```

3. Note the absolute path to the file `Harmonization.19JAN2017.GitHub.SourceFiles.R`

## How to Run

### From RStudio

Open R studio and the script `Harmonization.19JAN2017.GitHub.R`. Set these two variables:
```R
f.dir <- /phenotype/file/directory/ # from step 1 in Prereq's
source.file <- /my/local/git/repo/path/topmed-t2d-glycemia-public/methods/traitHarmonization/Harmonization.19JAN2017.GitHub.SourceFiles.R # from step 3 of Prereq's
```

*You may also need to comment out lines 15 and 16* 

### From the terminal

```bash
R --vanilla --args /phenotype/file/directory/ /my/local/git/repo/path/topmed-t2d-glycemia-public/methods/traitHarmonization/Harmonization.19JAN2017.GitHub.SourceFiles.R < Harmonization.19JAN2017.GitHub.R`
```
