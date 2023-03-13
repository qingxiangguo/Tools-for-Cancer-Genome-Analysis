# The usage of git command

```bash
# If you want to build a new repo
mkdir Allen_Guo_Rust_Learning_Note
cd Allen_Guo_Rust_Learning_Note 
git init
vim README.md

git clone https://github.com/qingxiangguo/qingxiangguo.git # clone a repo from git to local

git status

git add . # add all to the tree space

git commit -m "update"

git push origin main # The master has all been changed to main now

# The origin is the local place, the main(used to be master) is the main branch on the remote, you can 
# add and merge and switch between main branch and sub-branch

# If you want to delete something

rm -r ./*

git rm . # update the removal
```

If you init a local repo and want to push your local git repository to a remote repository on Github, you need to follow these steps:

Create a new repository on Github:

Go to https://github.com and log in to your account.
Click the "New" button on the top left corner and give your repository a name.
Click "Create repository" button.
Copy the remote repository URL:

After creating the repository, you will see a page with the repository details.
Copy the URL that starts with "https://github.com/username/repositoryname.git"
Add the remote repository to your local repository:

Open your terminal or command prompt and navigate to your local repository folder.
Type the following command to add the remote repository:

```bash
git remote add origin <remote_repository_URL> # Git creates a reference to the remote repository in your local repository. This reference is called a "remote" and it allows you to interact with the remote repository using Git commands.

git push origin main
```

