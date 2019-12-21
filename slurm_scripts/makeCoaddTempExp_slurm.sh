#!/bin/bash

# HEY YOU!!!
# BEFORE YOU SBATCH ANYTHING, YOU NEED TO MAKE A SKYMAP!!! e.g.,
# makeDiscreteSkyMap.py $repo --rerun $rerun --id object='Blind14A_04'^'Blind14A_09'^'Blind14A_10' --config skyMap.pixelScale=0.26

source /software/lsstsw/stack/loadLSST.bash
setup lsst_distrib

repo="/project/mrawls/hits2014"
rerun="2019_11:2019_12_coadds_full2"
configFile="/project/mrawls/ap_pipe/config/makeCoaddTempExp_goodSeeing.py"
moreConfig="select.nImagesMax=1000 select.maxPsfFwhm=5.33"
regularId="filter=g tract=0"
selectId="filter=g"  # object, i.e. field, is set by $1 (value in .conf file)

makeCoaddTempExp.py ${repo} --rerun ${rerun} -C ${configFile} --config ${moreConfig} --id ${regularId} --selectId ${selectId} $1 --clobber-versions
