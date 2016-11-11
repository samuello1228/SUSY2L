cd $1
dir=`pwd`
fileList=`ls`
allFiles=""

for file in $fileList
do 
	allFiles=$allFiles,$dir/$file
done

echo $allFiles