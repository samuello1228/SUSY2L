for FNAME in $( cat $1 ) ; do 
  root -b -q $ROOTCOREDIR/scripts/load_packages.C "NTUPPostProc.C+(\"$FNAME\")"
done
