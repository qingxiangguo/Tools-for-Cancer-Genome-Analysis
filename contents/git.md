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
