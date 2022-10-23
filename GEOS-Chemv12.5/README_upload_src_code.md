git rm --cached code_directory   # delete reference to submodule HEAD (no trailing slash)
git rm .gitmodules               # if you have more than one submodules,
			         # you need to edit this file instead of deleting!
rm -rf code_directory/.git       # make sure you have backup!!!!
git add code_directory/          # do not forget the slash!
git commit -m "remove submodule" # This will make the rm cached work!
git push origin main



**What if we want to add new files? Here goes useful commands.**
git status 
git add 
git commit -m "what has been changed?"
git push origin main
