#!/bin/bash

for i in `seq 0 9`
do
    mv mapp-"$i" mapp-0"$i"
    #mv mapp-"00$i" mapp-"$(($i+500))"
done

#for i in `seq 10 99`
#do
#    mv mapp-"$i" mapp-"$(($i+500))"
#done

