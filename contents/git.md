# Git Basics: A Quickstart Guide

## 1. Setting up a new repository

```bash
mkdir Allen_Guo_Rust_Learning_Note
cd Allen_Guo_Rust_Learning_Note 
git init
```

## 2. Cloning an existing repository

```bash
git clone https://github.com/qingxiangguo/qingxiangguo.git
```

## 3. Checking the status of your repository

```bash
git status
```

## 4. Making changes

Edit your files (e.g., using `vim README.md`), then:

```bash
git add .  # Add all changes to the staging area
git commit -m "Describe your changes here"
```

## 5. Pushing changes to a remote repository

```bash
git push origin main
```

## 6. Switching branches

To switch to an existing branch:

```bash
git checkout branch_name
```

To create and switch to a new branch based on your current branch:

```bash
git checkout -b new_branch_name
```

## 7. Fetching and merging changes from a remote repository

```bash
git fetch
git merge origin/main
```

Alternatively, you can use:

```bash
git pull origin main
```

## 8. Viewing commit history

To see the commit history of your current branch:

```bash
git log
```

To view the commits that are in `origin/main` but not in your local `main`:

```bash
git log main..origin/main
```

To see the commit history of the `main` branch in the `origin` remote:

```bash
git log origin/main
```

## 9. Deleting files and committing the removal

```bash
rm file_to_delete.txt
git rm file_to_delete.txt
git commit -m "Remove file_to_delete.txt"
```

## 10. Linking a local repository to a remote one

If you have initialized a local repository and want to connect it to a new remote repository:

1. Create a new repo on GitHub.
2. Copy the remote repository URL.
3. Add the remote repository to your local one:

```bash
git remote add origin remote_repository_URL
```

## 11. Pulling changes from a remote branch

```bash
git pull origin main
```

## Notes:

- `origin` is the default name given to the remote from which you have cloned your repository. It's essentially a shortcut to the remote repository's URL.
  
- `main` is often the default branch name (previously it was `master`). When you see something like `origin/main`, it refers to the `main` branch on the `origin` remote.

Below is the revised markdown guide, incorporating the error messages you provided and an explanation of the `git config pull.rebase false` setting, along with an overview of the differences between merge, rebase, and fast-forward:

# Handling Git Conflicts: A Merge Strategy Guide

Working with Git involves synchronizing your local changes with those in the remote repository. Conflicts arise when changes in both locations diverge. This guide will help you resolve conflicts using the merge strategy, which incorporates changes from the remote branch into your local branch.

## Initial Push Error

When you attempt to push your local changes, you might encounter the following error:

```
! [rejected]        main -> main (fetch first)
error: failed to push some refs to 'https://github.com/ylab-hi/SVoctopus.git'
hint: Updates were rejected because the remote contains work that you do not have locally.
```

This error indicates that there are commits in the remote repository that you don't have in your local branch.

## Attempting to Pull Changes

When you try to pull changes with `git pull origin main`, you may see an output like this:

```
remote: Enumerating objects: 15, done.
remote: Compressing objects: 100% (9/9), done.
From https://github.com/ylab-hi/SVoctopus
 * branch            main       -> FETCH_HEAD
fatal: Need to specify how to reconcile divergent branches.
```

Git is alerting you that your local branch has diverged from the remote branch, and you need to specify how to reconcile them.

## Configuring Git to Merge

To tell Git to use merging to reconcile divergent branches, you can set the configuration:

```bash
git config pull.rebase false
```

This command sets the `pull.rebase` option to `false`, which instructs Git to merge changes by default when you `pull`.

## Merge, Rebase, and Fast-Forward: What's the Difference?

- **Merge (`git config pull.rebase false`)**: Merging creates a new commit in your history that combines the changes from the remote branch with your local branch. This "merge commit" has two parents, one pointing to the last commit on your branch and one to the commits pulled from the remote.
  
- **Rebase (`git config pull.rebase true`)**: Rebasing rewrites your local branch's history by applying your changes on top of the remote branch's latest commit. This keeps a linear project history, as if you had made your local changes after the latest remote changes.

- **Fast-Forward (`git config pull.ff only`)**: A fast-forward merge can only occur when there have been no changes in your local branch that the remote branch does not already have. It simply moves your branch pointer to match the remote's, avoiding a merge commit.

## Pull and Merge Changes Again

After configuring with `git config pull.rebase false`, execute the pull command again:

```bash
git pull origin main
```

If there are conflicts, Git will pause and inform you which files need to be resolved manually.

## Resolve Conflicts and Commit

Open the conflicting files, resolve the changes, and then add and commit them:

```bash
git add .
git commit -m "Resolve merge conflicts"
```

## Push the Resolved Changes

With the conflicts resolved and the merge complete, you can now successfully push your changes:

```bash
git push origin main
```

## Give up all the local Changes and roll back to remote version
```bash
git checkout main
git fetch origin
git reset --hard origin/main
```
