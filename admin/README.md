# PRO-CF HOOMD-blue Admin

This repository was created to hold a modified version of HOOMD-blue v4.2.1 for colloid simulations in the PRO-CF group. This README details various admin functions that you may need to know about, including:

* [How to upgrade to a new HOOMD-blue version](/admin/README.md#how-to-upgrade-to-a-new-hoomd-blue-version)
* Tools for [comparing different HOOMD-blue versions](/admin/comparing-versions) 
* [How to add changes or make your own copy](/admin/README.md#how-to-make-changes) of this repository


## How to upgrade to a new HOOMD-blue version

Until we submit our changes to be integrated into the official HOOMD-blue release, you must manually copy all our modifications into any new release of HOOMD-blue that you wish to use.

It is important to do this periodically, so that we stay up-to-date

To do this, follow these steps:

1. Clone the latest stable release of HOOMD-blue.\n
**IMPORTANT**: this clone is linked to the original HOOMD-blue github repo and cannot be reuploaded to github as is.

2. Create a new directory called "unlinked-hoomd" and copy the hoomd-blue files into that directory as a new folder called "hoomd-v4.2.1" (this will allow you to always know which version you are using)

```bash
mkdir unlinked-hoomd
cp ../hoomd4.2.1-basic/hoomd-blue hoomd-v4.2.1
```

Now find the hidden ".git" folder and remove it

```bash
grep hoomd-v4.2.1 .git/*
rm -r [path-to-.git/].git/
```

3. Now go back to GitHub in your web browser and create a new repository caklled "hoomd4.2.1-mod" using the procf/general-template

4. Clone this new repo into local "hoomd4.2.1-mod" folder (you will now have a filepath `hoomd4.2.1-mod/hoomd4.2.1-mod/` but this is useful for keeping the virtual environment out of the github repo). Then cd into this repo and copy the "hoomd-v4.2.1" folder into it form the "unlinked-hoomd" directory called "hoomd-v4.2.1-props" (if the folder name ends in 4.2.1 it will be ignored by git)

```bash
mkdir hoomd4.2.1-mod
cd hoomd4.2.1-mod
git clone git@github.com:procf/hoomd4.2.1-mod.git
cd hoomd4.2.1-mod
cp [path-to-hoomd-v4.2.1] hoomd-v4.2.1-props
```

5. Push the new files to the Github. You should have normal folders containing all the hoomd-blue files for the version you want to modify! If any of the folders have a white arrow on them and cannot be openned... then you are still linked to the original hoomd-blue github and need to remove more .git files to delink it... sorry.

6. Now copy all your changes file-by-file from the old version to the latest version, paying attention to any larger level syntax changes that may have been implemented by the HOOMD-blue team.

7. In the new repo's main directory, update the README and add template scripts and other files.  

8. Test that the new version of hoomd-blue still runs colloidal sims correclty. You have now finally upgraded to the latest version!!
<br>
<br>
## How to Make Changes

You may need to make changes to this repository to add new features for everyone or to create an entirely new version for yourself.

1. [Jump to advice for adding changes for everyone](/README.md#how-to-add-changes)
2. Continue reading about making a copy

### How to Make a Copy

There are **two ways** to make your own, independent copy of this repositorywithout starting completely from scratch (as described above): 

1. Fork this repository! This will also copy the version history and retain links to this version. To create a fork: [Read about forks in the Github Docs](https://docs.github.com/en/get-started/quickstart/fork-a-repo)
2. Create a completely separate copy that is not linked to this repository. 

To create a copy that is not linked to this repository:

Download the basic version of HOOMD-blue as a tarball. Tarballs are available as GitHub release assets, as described on the [HOOMD-blue docs](https://hoomd-blue.readthedocs.io/en/stable/building.html#obtain-the-source). On the page for the version you want to install, copy the Download link and use it to curl the tar file.
```bash
curl -Lo hoomd-v4.2.1.tar.gz [link-to-tar, ex:https://github.com/glotzerlab/hoomd-blue/releases/download/v4.2.1/hoomd-v4.2.1.tar.gz]
tar -xvf hoomd-v4.2.1.tar.gz
rm hoomd-v4.2.1.tar.gz
```
If you are backing up to a new GitHub repository you may also need to remove the sphinx-doc folder
```bash
rm -r hoomd-v4.2.1/sphinx-doc
```
Delete the hoomd subfolder
```bash
rm -r hoomd-v4.2.1/hoomd
```
Use a copy of this repository's hoomd subfolder instead
```bash
cp -r [path-to-hoomd4.2.1-mod]/hoomd-v4.2.1-props/hoomd /hoomd-v4.2.1/hoomd
```
Rename the folder with the "-props" extension
```bash
mv hoomd-v4.2.1/ hoomd-v4.2.1-props
```
And, optionally, add back the sphinx-doc folder
```bash
cp -r [path-to-hoomd4.2.1-mod]/hoomd-v4.2.1-props/sphinx-doc /hoomd-v4.2.1-props/sphinx-doc
```
<br>

### How to Add Changes

To make changes and then merge them into this version:

1. Create a new branch to develop your changes
2. Label your final changes as described below
3. Ask someone else on the Team to review your code
4. Discuss with the Team to make sure the changes are thoroughly tested
5. Merge your finished branch into the main version, and update the changelog and full list of changed files

When adding mods to a file
- add a line in the header to note that the file has been changed at the current date (YYYY-MM)<br>
```C
///////~ PRO-CF in-house mods tagged //~ [PROCF-YYYY] ~///////
```
If this already exists, just update the date.

- At each change, give a brief explanation of your changes if you can, and then mark it with the tag referenced in the header
```C
added code /*commented out code*/ //~ optional explanation [PROCF-YYYY-MM]
```
Don't DELETE things, comment them out instead.

This will make the code easier to match with future versions of HOOMD-blue (when we need to transfer our mods over to a new veraion), and make the changes searchable for future users who may want to understand and/or further modify the modifications.
