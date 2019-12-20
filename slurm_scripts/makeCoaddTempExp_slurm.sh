#!/bin/bash

source /software/lsstsw/stack/loadLSST.bash
setup lsst_distrib

repo="/project/mrawls/hits2014"
rerun="2019_11:2019_12_coadds_full"
configFile="/project/mrawls/ap_pipe/config/makeCoaddTempExp_goodSeeing.py"
moreConfig="select.nImagesMax=1000 select.maxPsfFwhm=5.33"
regularId="filter=g tract=0"
selectId="filter=g"  # ccd is set by $1, i.e., value in .conf file

# BEFORE YOU RUN THIS, YOU'LL WANT TO MAKE A SKYMAP
# e.g., 
# makeDiscreteSkyMap.py ${repo} --rerun ${rerun} --id object='Blind14A_04'^'Blind14A_09'^'Blind14A_10' --config skyMap.pixelScale=0.26
makeCoaddTempExp.py ${repo} --rerun ${rerun} -C ${configFile} --config ${moreConfig} --id ${regularId} --selectId ${selectId} $1 --clobber-versions
