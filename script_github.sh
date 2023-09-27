#!/bin/bash

# Check the bullshit related to git
#git config --list

#user.name=thibauld
#user.email=thibauld_michel@protonmail.com
#color.ui=true
#credential.helper=cache
#core.repositoryformatversion=0
#core.filemode=true
#core.bare=false
#core.logallrefupdates=true
#remote.origin.url=https://github.com/ThibauldMichel/thesis_socotrana
#remote.origin.fetch=+refs/heads/*:refs/remotes/origin/*

# Start the shitshow
git init

# git add whatever it means
git add .

# git commit, no idea what it means
git commit -m "first commit" .

# Have a look at what is commited, whatever it means
#git log

# Add repo
#git remote add origin https://github.com/ThibauldMichel/Draft_PhD_thesis.git

# git push
#git push https://<GITHUB_ACCESS_TOKEN>@github.com/<GITHUB_USERNAME>/<REPOSITORY_NAME>.git

# TOKEN DE MERDE
#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
# git push https://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@github.com/ThibauldMichel/ROH-pipeline.git


#curl -H 'Authorization: token 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' https://github.com/ThibauldMichel/ROH-pipeline  


#curl -H 'Authorization: token <MYTOKEN>' ...

git remote remove origin
git remote add origin https://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@github.com/ThibauldMichel/ROH-pipeline.git
#git push
# git push --set-upstream origin main

git push --set-upstream origin main https://https://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@github.com/ThibauldMichel/ROH-pipeline.git

