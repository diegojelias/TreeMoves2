#!/bin/bash
mv git_newbranch.moves.sh ../

new_branch='oct11_imgay'

old_branch='oct9_ggm'

git checkout master
git branch -d $old_branch
git pull
git branch $new_branch
git checkout $new_branch
git push -u origin $new_branch

mv ../git_newbranch.moves.sh .
