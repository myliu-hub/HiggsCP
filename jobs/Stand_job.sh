#!/bin/bash

cd /scratchfs/bes/myliu/Liumy/HiggsCP/jobs
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

InputDir_pre=/cefs/data/stdhep/CEPC240/higgs
InputDir=dirflag
InputFile=fileflag

OutputDir_pre=/cefs/higgs/myliu/Higgs_CP_Data/higgs
OutputDir=dirflag
OutputFile=${InputFile}.root
/scratchfs/bes/myliu/Delphes/DelphesSTDHEP /scratchfs/bes/myliu/Delphes/cards/delphes_card_CEPC.tcl $OutputDir_pre/$OutputDir/$OutputFile $InputDir_pre/$InputDir/$InputFile
