#!/bin/bash
AllFileList=/scratchfs/bes/myliu/Liumy/HiggsCP/jobs/txt/File_higgs.txt
echo $AllFileList
for FileDir in `cat $AllFileList`
do
    FileList=$FileDir
    echo $FileList
    Dir=`basename $FileList`
    OutputDir_pre=/cefs/higgs/myliu/Higgs_CP_Data/higgs
    Shell=/cefs/higgs/myliu/Higgs_CP_Data/Log

    [ -d $OutputDir_pre/$Dir ] || mkdir $OutputDir_pre/$Dir
    [ -d $Shell/$Dir ] || mkdir $Shell/$Dir
    for file in `cat $FileList`
    do
        name=`basename $file`
        echo $Dir $name
        sed -e 's/dirflag/'$Dir'/g' -e 's/fileflag/'$name'/g' Stand_job.sh > $Shell/$Dir/${name}.sh
        chmod +x  $Shell/$Dir/${name}.sh
        hep_sub -g physics $Shell/$Dir/${name}.sh -e $Shell/$Dir/ -o $Shell/$Dir/
    done
done
