#!/bin/bash
param_1=$1
param_2=$2

for i in $(ls $param_1);do
    fichier=`echo $param_1"/"$i` 
    prefix=`echo $param_2"/"$i`
    split -b 50MB $fichier $prefix
done
