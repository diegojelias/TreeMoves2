#!/bin/bash

f="Cluster_maketrees_25_1000.py"
base=`basename $f .py`
for x in {1..10}
do
echo $base"_"$x
cp $f "temp.py"
sed -i.tmp "s/taco/$x/g" "temp.py"
python temp.py
rm "temp.py"
rm "temp.py.tmp"
done 

f="Cluster_maketrees_50_100.py"
base=`basename $f .py`
for x in {1..10}
do
echo $base"_"$x
cp $f "temp.py"
sed -i.tmp "s/taco/$x/g" "temp.py"
python temp.py
rm "temp.py"
rm "temp.py.tmp"
done 

f="Cluster_maketrees_50_1000.py"
base=`basename $f .py`
for x in {1..10}
do
echo $base"_"$x
cp $f "temp.py"
sed -i.tmp "s/taco/$x/g" "temp.py"
python temp.py
rm "temp.py"
rm "temp.py.tmp"
done 

f="Cluster_maketrees_25_100.py"

base=`basename $f .py`
for x in {1..10}
do
echo $base"_"$x
cp $f "temp.py"
sed -i.tmp "s/taco/$x/g" "temp.py"
python temp.py
rm "temp.py"
rm "temp.py.tmp"
done 

