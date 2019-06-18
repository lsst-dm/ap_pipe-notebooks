Ingestion of HiTS instcals
==========================

Memo for creating repo on lsst-dev for DM-19525 instcal processing. 
We performed 2 processing sessions of the HiTS 2015 instcals using a calexp 2014 template.
Locations:
::

    1: /project/gkovacs/hits_proc2019-06-20 (earlier: hits_instcals)
    2: /project/gkovacs/hits_proc2019-07-03

The 1st session processed 2 fields with 2+1 visits:

Source files:

- Linked to Eric's download: /project/gkovacs/hits2014; /project/gkovacs/hits2015
- /project/gkovacs/hits_cpline_dl

NOAO archive downloads can be best filtered by the Proposal IDs: 2015A-0608 and 2014A-0608

Processing session 2019-06-18
=============================

Following instructions on obs_decam doc pages how to prepare repo for
``processCcd``.  Separate files into ``instcal``, ``dqmasks``,
``wtmaps`` directories.
::

    mkdir repo
    echo lsst.obs.decam.DecamMapper > repo/_mapper

Use ``ingestImagesDecam.py`` to ingest data files. It automatically
replaces all instances of the substring ``instcal`` in the given full absolute path 
to look for the dqmasks and wtmaps. Do not use ``instcal`` anywhere else in the full path.
Repeat command for ``hits2014`` and ``hits2015`` directories.
::

    ingestImagesDecam.py ../repo/ --filetype instcal --mode=link instcals/*.fits.fz

Create rerun folder and initalize association db:
::

    mkdir -p rerun/ap_pipe_2019-06-13
    make_ppdb.py -c ppdb.db_url=sqlite:///repo/rerun/ap_pipe_2019-06-13/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"

Copy of reference catalogs:
::

    mkdir ref_cats
    mkdir gaia
    mkdir pan-starrs
    cd pan-starrs/
    tar -xvzf /project/mrawls/ap_verify_ci_hits2015/refcats/ps1_HiTS_2015.tar.gz
    cd ../gaia/
    tar -xvzf /project/mrawls/ap_verify_ci_hits2015/refcats/gaia_HiTS_2015.tar.gz 

ap_pipe notes
-------------

ap_pipe command line arguments have hard wired datatype ``raw`` for the ``--id`` dataset field. This means that ``runDataRef`` is not called if the repo does not contain any ``raw`` data. This could be improved by using DynamicDataset in the argument parser. ``processCcd` is not affected by this because it uses the dataset type from its ``isr`` subtask and it is redirected to DecamNullIsr for instcals.

Otherwise ap_pipe supports reusing of processCcd output and picking up ``calexp``-s, though as of 2019-06-25 this works only for writeable
calexp locations (likely trivial bug).

Running of processCcd
---------------------
::

    [gkovacs@lsst-dev03 hits_instcals]$ processCcd.py repo/ --rerun proc_2019-06-18 -C config/ processCcd_2019-06-18.py --id |& tee -a procCcd_2019-06-18.log

Running of imageDifference
--------------------------
::

    [gkovacs@lsst-dev03 hits_instcals]$ bash run_imageDifference_2019-06-20.sh |& tee -a imageDifference_2019-06-20.log

Running for each visit with its matching template image:
::

    imageDifference.py repo/ --rerun proc_2019-06-18:imgDiff_2019-06-20 -C config/imageDifference_2019-06-20.py --id filter=g visit=411420^419802 --templateId visit=289450
    imageDifference.py repo/ --rerun proc_2019-06-18:imgDiff_2019-06-20 -C config/imageDifference_2019-06-20.py --id filter=g visit=411371 --templateId visit=289449

Running of ap_association
--------------------------
::

    [gkovacs@lsst-dev03 hits_instcals]$ make_ppdb.py -c ppdb.db_url=sqlite:///repo/rerun/imgDiff_2019-06-20/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"
    [gkovacs@lsst-dev03 hits_instcals]$ python3 ../assoc_wrapper/run_assoc_2019-06-21.py |& tee -a run_assoc_2019-06-26.log
