#!/bin/bash
for i in $(seq -f "%02.f" 0 99)
do 
cd mapp-${i}
 cat fort.200 | awk '{print $2}' | awk -v var="$i" '{if( $1 == "Infinity") print NR, $1 , var}'
cd ../
done

