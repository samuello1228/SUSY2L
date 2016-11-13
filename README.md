*Basic instructions about how to submit updates*

*1. Create a branch from the master*

```sh
git clone ssh://git@gitlab.cern.ch:7999/hku/SUSY2L.git # takes the master

git checkout -b "the name of your branch" # switches locally to a branch

git push -u origin "the name of your branch" # pushed the branch on GitLab
```

*2. Develop on your branch*

(Here you can do whatever you want and updated the package as you wish)
```sh
git add example_of_file # this is a file where you implemented some changes

git commit -m "commit message"

git status # (you have to check you are on your branch)

git push # this puts the commits on the web in your branch only (the master is unaffected)
```
*3. Make your tests*

*4. Request a merge into the master*

Assign it to somebody that will check your changes and make independent tests

(Done via web interface or terminal commands)

*5. Then you can tag this master*

(Done via web interface or terminal commands)