#! /bin/bash


i=1
list=`ls *.xyz`
for item in $list
do
 tail -n +3 $item > $item.mid
 cat sp_header $item.mid sp_footer >> $item.inp
 rm -f $item.mid
 ID=`printf "%0*d\n" 3 $i`
 mv $item.inp "q$ID.$item.inp"
 echo "file $i: $item --> q$ID.$item"
 let i=i+1
done


