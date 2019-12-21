#!/bin/bash

source /software/lsstsw/stack/loadLSST.bash
setup lsst_distrib

repo="/project/mrawls/hits2014"
rerun="2019_11:2019_12_coadds_full2"
assembleConfig="doInterp=True doNImage=True"
assembleId="filter=g tract=0"  # patch is set by $1, i.e., value in .conf file
selectId="filter=g"

assembleCoadd.py ${repo} --rerun ${rerun} --warpCompareCoadd --selectId ${selectId} --id ${assembleId} $1 --config ${assembleConfig} --clobber-versions --clobber-config
