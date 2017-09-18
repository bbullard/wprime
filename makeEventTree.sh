#!/bin/bash/

export WORKDIR="/export/home/bbullard/thesis"
MCDIR="/export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/mc"
DATADIR="/export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/data"
ls -d -1 ${MCDIR}/*Wplus*Ntuple*/ > ${WORKDIR}/W.txt
ls -d -1 ${MCDIR}/*Wmin*Ntuple*/ >> ${WORKDIR}/W.txt
ls -d -1 ${MCDIR}/*ttbar*Ntuple*/ > ${WORKDIR}/Top.txt
ls -d -1 ${MCDIR}/*top*Ntuple*/ >> ${WORKDIR}/Top.txt
ls -d -1 ${MCDIR}/*DY*Ntuple*/ > ${WORKDIR}/DY.txt
ls -d -1 ${MCDIR}/*Zmumu*Ntuple*/ >> ${WORKDIR}/DY.txt
ls -d -1 ${MCDIR}/*Ztautau*Ntuple*/ >> ${WORKDIR}/DY.txt
ls -d -1 ${MCDIR}/*Sherpa*Ntuple*/ > ${WORKDIR}/Dib.txt
ls ${DATADIR}/user.dhaliwal.Wprime.v9.data16_13TeV.002*Ntuple*/* > ${WORKDIR}/data.txt
ls ${DATADIR}/user.dhaliwal.Wprime.v9.data16_13TeV.003*Ntuple*/* >> ${WORKDIR}/data.txt
ls ${DATADIR}/user.dhaliwal.Wprime.v9.data15*Ntuple*/* >> ${WORKDIR}/data.txt
ls -d -1 /export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/signal/*2000_Ntuple*/ > ${WORKDIR}/signal2.txt
ls -d -1 /export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/signal/*3000_Ntuple*/ > ${WORKDIR}/signal3.txt
ls -d -1 /export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/signal/*4000_Ntuple*/ > ${WORKDIR}/signal4.txt
ls -d -1 /export/home/dhaliwal/Files/W-prime/ntuples/v9_flxplus/signal/*5000_Ntuple*/ > ${WORKDIR}/signal5.txt
#sed -i '/301154/d' ${WORKDIR}/W.txt
#sed -i '/301134/d' ${WORKDIR}/W.txt

#cp ${WORKDIR}/W.txt ${WORKDIR}/Files/ntuples/selection/Wcutflow.txt 
#cp ${WORKDIR}/Top.txt ${WORKDIR}/Files/ntuples/selection/Topcutflow.txt 
#cp ${WORKDIR}/DY.txt ${WORKDIR}/Files/ntuples/selection/DYcutflow.txt
#cp ${WORKDIR}/Dib.txt ${WORKDIR}/Files/ntuples/selection/Dibcutflow.txt 
#cp ${WORKDIR}/data.txt ${WORKDIR}/Files/ntuples/selection/datacutflow.txt

#cp ${WORKDIR}/W.txt ${WORKDIR}/Files/ntuples/selection/WQCDcutflow.txt 
#cp ${WORKDIR}/Top.txt ${WORKDIR}/Files/ntuples/selection/TopQCDcutflow.txt 
#cp ${WORKDIR}/DY.txt ${WORKDIR}/Files/ntuples/selection/DYQCDcutflow.txt
#cp ${WORKDIR}/Dib.txt ${WORKDIR}/Files/ntuples/selection/DibQCDcutflow.txt 

export TXTLOC="${WORKDIR}/W.txt"
export RFILE="WFullSelection.root"
export MODE="W"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/Top.txt"
export RFILE="TopFullSelection.root"
export MODE="Top"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/DY.txt"
export RFILE="DYFullSelection.root"
export MODE="DY"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/Dib.txt"
export RFILE="DibFullSelection.root"
export MODE="Dib"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/data.txt"
export RFILE="dataFullSelection.root"
export MODE="data"
#condor_submit ${WORKDIR}/makeEventTree.sub


#Signal MC

export TXTLOC="${WORKDIR}/signal3.txt"
export RFILE="signal3FullSelection.root"
export MODE="signal3"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/signal2.txt"
export RFILE="signal2FullSelection.root"
export MODE="signal2"
#condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/signal4.txt"
export RFILE="signal4FullSelection_LB.root"
export MODE="signal4_LB"
condor_submit ${WORKDIR}/makeEventTree.sub

export TXTLOC="${WORKDIR}/signal5.txt"
export RFILE="signal5FullSelection.root"
export MODE="signal5"
#condor_submit ${WORKDIR}/makeEventTree.sub


#QCD Background estimation
export TXTLOC="${WORKDIR}/W.txt"
export RFILE="WQCD.root"
export MODE="WQCD"
#condor_submit ${WORKDIR}/makeEventTreeQCD.sub

export TXTLOC="${WORKDIR}/Top.txt"
export RFILE="TopQCD.root"
export MODE="TopQCD"
#condor_submit ${WORKDIR}/makeEventTreeQCD.sub

export TXTLOC="${WORKDIR}/DY.txt"
export RFILE="DYQCD.root"
export MODE="DYQCD"
#condor_submit ${WORKDIR}/makeEventTreeQCD.sub

export TXTLOC="${WORKDIR}/Dib.txt"
export RFILE="DibQCD.root"
export MODE="DibQCD"
#condor_submit ${WORKDIR}/makeEventTreeQCD.sub

export TXTLOC="${WORKDIR}/data.txt"
export RFILE="dataQCD.root"
export MODE="dataQCD"
#condor_submit ${WORKDIR}/makeEventTreeQCD.sub

