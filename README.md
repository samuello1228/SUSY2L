## Basic instructions about how to submit updates

### 0. Configure git on your computer
```sh
git config --global user.name "Your name"
git config --global user.email "Your email"
git config --global core.editor vi
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.up rebase
git config --global alias.ci commit
```

### 1. Create a branch from the master

```sh
# Clone the repository to your computer
git clone ssh://git@gitlab.cern.ch:7999/hku/SUSY2L.git # via ssh
git clone https://:@gitlab.cern.ch:8443/hku/SUSY2L.git/ # via Kerberos # do kinit user@CERN.CH first

# Make a new branch and switch to it
git checkout -b "the name of your branch"

git push -u origin "the name of your branch" # pushed the branch on GitLab
```

### 2. Develop on your branch

(Here you can do whatever you want and updated the package as you wish)
```sh
git add example_of_file # this is a file where you implemented some changes

git commit -m "commit message"

git status # (you have to check you are on your branch)

git push # this puts the commits on the web in your branch only (the master is unaffected)
```

### 3. Get a existing branch from the origin

```sh
git checkout -b newBranch origin/oldBranch
```

### 4. Request a merge into the master

Assign it to somebody that will check your changes and make independent tests

(Done via [web interface](https://gitlab.cern.ch/hku/SUSY2L/branches) or terminal commands)

### 5. Then you can tag this master

(Done via web interface or terminal commands)
