#!/bin/bash
> DataSample.txt
name=physics_Main
zero=00

exec 3< DataLumi.txt
while read a b <&3
do
  c=${zero}$a
  echo $c ${name} $b>> DataSample.txt
done

cat DataSample.txt
