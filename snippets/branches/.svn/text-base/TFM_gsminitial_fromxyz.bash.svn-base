#/bin/bash

# take all xyz files, make q-chem inputs

i=1
ID=`printf "%0*d\n" 3 $i`
list1=`ls 3?.xyz`
echo $list1
for item in $list1
do
  var1="${item%????}" 
  list2=`ls $var1-???.xyz`
  for item2 in $list2
    do
    var2="${item2%????}"
    ID=`printf "%0*d\n" 3 $i`
    cat $item $item2 >> scratch/q$ID.$var1.TS.$var2.xyz
    list3=`ls $var2-???.xyz`
    let i=i+1
    for item3 in $list3
      do
      var3="${item3%????}"
      ID=`printf "%0*d\n" 3 $i`
      cat $item2 $item3 >> scratch/q$ID.$var2.TS.$var3.xyz
      list4=`ls $var3-???.xyz`
      let i=i+1
      for item4 in $list4
	do
	var4="${item4%????}"
        ID=`printf "%0*d\n" 3 $i`
        cat $item3 $item4 >> scratch/q$ID.$var3.TS.$var4.xyz
	list5=`ls $var4-???.xyz`
	let i=i+1
        for item5 in $list5
	  do
	  var5="${item5%????}"
          ID=`printf "%0*d\n" 3 $i`
	  cat $item3 $item4 >> scratch/q$ID.$var3.TS.$var4.xyz
	  list5=`ls $var5-???.xyz`
	  let i=i+1
 	done
      done
    done      
  done
done  


#  ID=`printf "%0*d\n" 3 $i`
#  echo "file $i: $item --> q$ID.$item.inp"
#  csplit $item.xyz 3 56 
#  cat header xx01 footer > scratch/q$ID.$item.inp 
#
#  let i=i+1
#done
