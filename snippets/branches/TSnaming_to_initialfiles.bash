#! /bin/bash

listxyz=`ls *.xyz`
i=1
for item in $listxyz
do
  ID=`printf "%0*d\n" 4 $i`
#  mv $item "q$ID.$item"
#  mv $item "initial$ID.xyz"
  echo "file $i: $item --> initial$ID.xyz" >> tsqfilenamechange.txt
  let i=i+1
done
