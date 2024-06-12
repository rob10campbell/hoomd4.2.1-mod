# Adding New Changes 

You may need to make changes to this repository to add new features for everyone or to create an entirely new version for yourself.

Here are some tips and best-practices for integrating new changes into hoomd4.2.1-mod:

1. [Jump to advice for adding changes for everyone](/README.md#how-to-add-changes)
2. Continue reading about making a personal copy

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

To make changes and then merge them into this version, follow a modified version of the standard Github Flow ([more details here](https://docs.github.com/en/get-started/using-github/github-flow#following-github-flow))

1. Create a new branch to develop your changes
2. Label your final changes as described below
3. Merge with the main branch to make sure your changes work with the latest version
4. Create a [detailed pull request](/admin/pull-request-template.md) so someone else on the Team can review your code
5. Discuss with the Team to make sure the changes are thoroughly tested
6. Merge your finished branch into the main version, update the changelog and full list of changed files, and delete your branch

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
