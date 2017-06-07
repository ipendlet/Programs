#/bin/bash


list1=`ls *.out`
for item in $list1
  do
  linenumber=`grep -n "END OF MOLDEN-FORMATTED" $item | cut -d: -f1`
  let firstline=$linenumber-56
  let lastline=$linenumber-1
  csplit $item $firstline $lastline
  mv xx01 "${item:5:${#item}-13}"
done
