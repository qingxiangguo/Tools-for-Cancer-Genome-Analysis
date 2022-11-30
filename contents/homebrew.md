# The installation and usage of Homebrew
# 1. About
Homebrew is the second mac APP store. Like conda and mamba.
# 2. Installation

Paste that in a macOS Terminal or Linux shell prompt.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Add to environmental variant
```
echo '# Set PATH, MANPATH, etc., for Homebrew.' >> /Users/qgn1237/.profile
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/qgn1237/.profile
eval "$(/opt/homebrew/bin/brew shellenv)"
```

Make sure  you have the newest Homebrew

```
brew upgrade
```
# 3. Usage

Install software with GUI, you should add --cask, or it will only be in command line.
```
brew install --cask iterm2 # local terminal

brew install --cask zotero

brew install --cask rectangle # for screen management
```




