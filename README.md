# NOTES ON TRAIT DATA

This is a *public* repositiory and as such, trait data should *never* be uploaded. This repository is only intended to hold scripts for trait harmonization, not actual trait files or any data that can be construed as violating consent/permission/privacy aggreements set through TOPMed.

Certain restrictions are in place to not allow specific file types to be uploaded. However, please make sure that not trait data is included in the comments of any script.

# Instructions for contributing 

Whenever a user would like to contribute (edit or add files), it should be done on a new branch. Assuming the repository has been cloned locally, follow these steps.

## Creating a new branch

```bash

git checkout -b my_branch

```

## Modifying code

Once a new branch has been created, make any code changes to the files in the local repository.

## Staging changes to commit

```bash
git status # this will show what files have been changed
```

```bash
git add my_changed_file.R # this stage the changes to be committed, only the changes that you add will be pushed back to the remote repository
```

```bash
git commit -m "a message describing my changes" # this will add a helpful message to the changes
```

## Pushing changes to the remote repository

```bash
git push origin my_branch # this will push the changes to the remote repository
```

## Pull request and merging changes into the master branch

Once the changes have been pushed, navigate to the [repository](https://github.com/manning-lab/topmed-t2d-wg-trait-harmonization/pulls) in a web browser. Click *New pull request* and select *my_branch* to be merged into *master*. Then click create pull request. On the next page, click *Merge*. This will move your changes into the *master* branch.

