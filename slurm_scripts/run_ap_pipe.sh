#!/bin/bash

# This script is called by run_ap_pipe.conf, which is called by run_ap_pipe.sl

# HEY YOU!!!
# BEFORE YOU SBATCH ANYTHING, YOU NEED TO MAKE AN EMPTY APDB!!! e.g.,
# make_apdb.py -c apdb.isolation_level=READ_UNCOMMITTED -c apdb.db_url="sqlite:////project/mrawls/hits2015/rerun/$rerun/association.db"

# ALSO!!!
# BEFORE YOU SBATCH ANYTHING, YOU MIGHT WANT TO LOOK AT THE CONFIG SETTINGS!!!
# the optional custom config file is called ap_pipe_config.py

# Set up the stack and necessary ap_packages
source /software/lsstsw/stack/loadLSST.bash
setup lsst_distrib

# Config stuff
rerun="cw_2019_12"
calib="/project/mrawls/hits2015/calib5"
template="/project/mrawls/hits2014/rerun/2019_12_coadds_full"
#template="/project/sullivan/hits2014B/cwcoadds_processed2" # high-res
#template="/project/mrawls/hits2014/coadds_processed2"  # old standard
db_config=(-c apdb.db_url="sqlite:////project/mrawls/hits2015/rerun/$rerun/association.db" -c apdb.isolation_level="READ_UNCOMMITTED" -c apdb.connection_timeout=240)
#more_config=(-C /project/mrawls/hits2015/ap_pipe_config.py)
more_config=(-c ccdProcessor.calibrate.photoCal.match.referenceSelection.magLimit.fluxField="i_flux" -c ccdProcessor.calibrate.photoCal.match.referenceSelection.magLimit.maximum=22.0)

# Command to run
ap_pipe.py /project/mrawls/hits2015 --calib ${calib} --template ${template} --rerun ${rerun} "${db_config[@]}" "${more_config[@]}" --longlog --id $*
