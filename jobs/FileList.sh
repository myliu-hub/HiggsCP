#!/bin/bash

Dirlist=txt/Dir_higgs.txt
Dir_pre=/cefs/data/stdhep/CEPC240/higgs/
for dir in `cat $Dirlist`
do
   ls $Dir_pre/$dir/* >./Higgs/$dir
done
