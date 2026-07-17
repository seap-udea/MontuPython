# Agent Rules

## make stash
When I (the user) request to "make stash" or ask you to pull changes from a remote repository without losing my local changes, follow this exact procedure:
1. Run `git stash` to save the local changes safely.
2. Run `git pull` to fetch and integrate the latest remote changes.
3. Run `git stash pop` to re-apply the local changes on top of the updated branch.
4. If the instruction also includes pushing the changes (e.g. "y después haz push"), proceed to stage (`git add`), commit (`git commit -m "..."`), and push (`git push`) the changes.
