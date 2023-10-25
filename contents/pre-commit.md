# Pre-commit Tutorial

`pre-commit` is a tool used to manage and maintain pre-commit hooks. These hooks help ensure that your code meets certain standards before it's committed to a repository. Here's how to set up and use `pre-commit`.

## **Step 1: Inside a Git Repository**
Ensure you are within an initialized Git repository. If you're not, you can initialize one using:
```bash
git init
```

## **Step 2: Create a `.pre-commit-config.yaml` File**
This file defines the hooks you'd like to use. You can either manually create this file or use the `pre-commit sample-config` command to generate a sample configuration. Then, adjust it according to your needs.
```bash
pre-commit sample-config > .pre-commit-config.yaml
```
Within this file, you can specify multiple sources for hooks and the particular hooks you wish to use.

## **Step 3: Install the `pre-commit` Hooks**
By executing the following command, `pre-commit` will be installed into the `.git/hooks/` directory, ensuring it runs every time an attempt to commit is made.
```bash
pre-commit install
```

## **Step 4: Manually Running Hooks**
Before making a commit, you might want to manually run the hooks to ensure everything is functioning as expected. This is especially useful when adding new hooks or making significant changes.
```bash
pre-commit run --all-files
```

## **Step 5: Making a Git Commit**
Now, every time you try to make a Git commit, `pre-commit` will automatically run the hooks you've specified in `.pre-commit-config.yaml`. 
```bash
git add <your-changes>
git commit -m "Your commit message"
```
If any of the hooks detect issues, they will halt the commit and report the problems. You can then fix these issues and try to commit again.

---

## **Conclusion**
In essence, `pre-commit` is a powerful tool that ensures code meets certain standards before it's committed to a repository. By automating these checks, it assists teams in maintaining code quality and catching potential issues early on.