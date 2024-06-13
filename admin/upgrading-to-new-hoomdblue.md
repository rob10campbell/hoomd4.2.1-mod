# Upgrading to a New HOOMD-blue Version

This repository was created to hold a modified version of HOOMD-blue v4.2.1 for colloid simulations in the PRO-CF group. In the future you may want to integrate these changes into a newer version of HOOMD-blue; however, this is complicated and time consuming. HOOMD-blue is constantly improving, but in the past they have also removed features that were essential to our work. Be aware that our changes are specific to our research, and may not easily integrate with other updates.

Jump to:

* [How to upgrade to a new HOOMD-blue version](/admin/upgrading-to-new-hoomdblue.md#how-to-upgrade-to-a-new-hoomd-blue-version)
* Tools for [comparing different HOOMD-blue versions](/admin/comparing-versions) 


## How to upgrade to a new HOOMD-blue version

Until we submit our changes to be integrated into the official HOOMD-blue release, you must manually copy all our modifications into any new release of HOOMD-blue that you wish to use.

It is important to do this periodically, so that we stay up-to-date

To do this:

1. Go to GitHub in your web browser and create a new repository for your version (for example "hoomd4.2.1-mod") using the procf/general-template

2. On your local computer, create a new folder for your version of the software: "hoomd4.2.1-mod" and then cd into it and clone the git repo you just made.<br>
(you will now have a filepath `hoomd4.2.1-mod/hoomd4.2.1-mod/` but this is useful for keeping the virtual environment out of the github rep)

3. Copy all the administrative documents from the previous version into this new git repo and push the changes

4. Now download the latest version of HOOMD-blue from the tarball link provided on the [HOOMD-blue website](https://hoomd-blue.readthedocs.io/en/v4.2.1/building.html#obtain-the-source) <br>
**IMPORTANT**: If you clone the GitHub repo instead, then the new copy of HOOMD-blue will be linked to the original HOOMD-blue github repo and cannot be reuploaded to GitHub as a separate set of files.

```bash
# replace hoomd-v4.2.1 with the version you are installing
curl -Lo hoomd-v4.2.1.tar.gz https://github.com/glotzerlab/hoomd-blue/releases/download/v4.2.1/hoomd-4.2.1.tar.gz
tar -xvf hoomd-v4.2.1.tar.gz   
```

5. Rename this file to add the "-procf" extension so it will never be ignored by the .gitignore file AND you always know what version you last updated (e.g. hoomd-v4.2.1 -> hoomd-v4.2.1-procf)

5. You can now copy our changes into the new version, make any additional modifications needed to keep them working, and update all the install/update and template scripts accordingly!
<br>
<br>
If for some reason you want to clone the git repo, or the tarball is still linked to the original github, you can break the link by removing the .git folder:

4. Create a new directory called "unlinked-hoomd" and copy the hoomd-blue files into that directory as a new folder called "hoomd-v4.2.1" (this will allow you to always know which version you are using)

```bash
mkdir unlinked-hoomd
cp ../hoomd4.2.1-basic/hoomd-blue hoomd-v4.2.1
```

Now find the hidden ".git" folder and remove it

```bash
grep hoomd-v4.2.1 .git/*
rm -r [path-to-.git/].git/
```

5. Copy the "hoomd-v4.2.1" folder from the "unlinked-hoomd" directory into your new git repo direcotry, and rename it "hoomd-v4.2.1-procf" (if the folder name ends in 4.2.1 it will be ignored by git)

6. Push the new files to the Github. You should have normal folders containing all the hoomd-blue files for the version you want to modify! If any of the folders have a white arrow on them and cannot be openned... then you are still linked to the original hoomd-blue github and need to remove more .git files to delink it... sorry.

7. Now copy all your changes file-by-file from the old version to the latest version, paying attention to any larger level syntax changes that may have been implemented by the HOOMD-blue team.
