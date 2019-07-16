Processing session 2019-07-03
=============================

All helper scripts are stored in
``ap_pipe-notebooks/false_positives/DM-19525_HiTS_instcalProc_calexpTemplate``.

Create the folders
------------------
::

    [gkovacs@lsst-dev03 Blind14A_10]$ mkdir instcals
    [gkovacs@lsst-dev03 Blind14A_10]$ mkdir wtmaps
    [gkovacs@lsst-dev03 Blind14A_10]$ mkdir dqmasks

Find the 2014 visit we want to use
------------------------------------
::

    [gkovacs@lsst-dev03 Blind14A_10]$ python3 print_visit_nums_fits.py *fits* | grep 289450
        289450   tu2203328.fits.fz
        289450   tu2203517.fits.fz
        289450   tu2203700.fits.fz

Distribute input instcals into dirs for ingestion
--------------------------------------------------

Blind14A_10
~~~~~~~~~~~~
::

    [gkovacs@lsst-dev03 Blind14A_10]$ python3 dist_image_weightmap_dqmap.py tu2203328.fits.fz tu2203517.fits.fz tu2203700.fits.fz
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_10/tu2203328.fits.fz ../wtmaps/
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_10/tu2203517.fits.fz ../instcals/
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_10/tu2203700.fits.fz ../dqmasks/
    [gkovacs@lsst-dev03 Blind14A_10]$ python3 dist_image_weightmap_dqmap.py tu2203328.fits.fz tu2203517.fits.fz tu2203700.fits.fz | bash

Blind14A_09
~~~~~~~~~~~~
::

    [gkovacs@lsst-dev03 Blind14A_09]$ python3 print_visit_nums_fits.py *fits* | grep 289449
        289449   c4d_140303_050532_ood_g_a1.fits.fz
        289449   c4d_140303_050532_ooi_g_a1.fits.fz
        289449   c4d_140303_050532_oow_g_a1.fits.fz
        289449   tu2202825.fits.fz
        289449   tu2203670.fits.fz
        289449   tu2203702.fits.fz
    [gkovacs@lsst-dev03 Blind14A_09]$ python3 dist_image_weightmap_dqmap.py tu2202825.fits.fz tu2203670.fits.fz  tu2203702.fits.fz
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_09/tu2202825.fits.fz ../instcals/
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_09/tu2203670.fits.fz ../dqmasks/
    ln -s /project/gkovacs/hits_instcal_dl/Blind14A_09/tu2203702.fits.fz ../wtmaps/
    [gkovacs@lsst-dev03 Blind14A_09]$ python3 dist_image_weightmap_dqmap.py tu2202825.fits.fz tu2203670.fits.fz   tu2203702.fits.fz | bash

Blind14A_04
~~~~~~~~~~~~
::

    [gkovacs@lsst-dev03 Blind14A_04]$ python3 print_visit_nums_fits.py *fits* | grep 289444
        289444   c4d_140303_044946_ood_g_a1.fits.fz
        289444   c4d_140303_044946_ooi_g_a1.fits.fz
        289444   c4d_140303_044946_oow_g_a1.fits.fz
        289444   tu2202855.fits.fz
        289444   tu2202994.fits.fz
        289444   tu2203813.fits.fz
    [gkovacs@lsst-dev03 Blind14A_04]$  python3 dist_image_weightmap_dqmap.py tu2202855.fits.fz tu2202994.fits.fz tu2203813.fits.fz | bash


Blind15A_26 Blind15A_40 Blind15A_42
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    [gkovacs@lsst-dev03 hits_cpline_dl]$ cd Blind15A_26/
    [gkovacs@lsst-dev03 Blind15A_26]$ python3 dist_image_weightmap_dqmap.py *g_v1.fits.fz | bash
    [gkovacs@lsst-dev03 Blind15A_26]$ cd ../Blind15A_42/
    [gkovacs@lsst-dev03 Blind15A_42]$ python3 dist_image_weightmap_dqmap.py *g_v1.fits.fz | bash
    [gkovacs@lsst-dev03 Blind15A_42]$ cd ../Blind15A_40/
    [gkovacs@lsst-dev03 Blind15A_40]$ python3 dist_image_weightmap_dqmap.py *g_v1.fits.fz | bash
    [gkovacs@lsst-dev03 Blind15A_40]$

Check filters for old name style 2014A (template) files
----------------------------------------------------------
::

    [gkovacs@lsst-dev03 instcals]$ fitsheader -k FILTER -e 0 tu2202825.fits.fz tu2202855.fits.fz tu2203517.fits.fz
        # HDU 0 in tu2202825.fits.fz:
        FILTER  = 'g DECam SDSS c0001 4720.0 1520.0' / Unique filter identifier
        # HDU 0 in tu2202855.fits.fz:
        FILTER  = 'g DECam SDSS c0001 4720.0 1520.0' / Unique filter identifier
        # HDU 0 in tu2203517.fits.fz:
        FILTER  = 'g DECam SDSS c0001 4720.0 1520.0' / Unique filter identifier

Ingestion
---------

Ingestion input abspath MUST NOT contain "instcal", "wtmap" "dqmask" at other locations.
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ ingestImagesDecam.py repo/ --filetype instcal --mode=link ../hits_cpline_dl/instcals/*.fits.fz |& tee -a ingest_2019-07-03.log


Reference catalogs
------------------
::

    [gkovacs@lsst-dev03 repo]$ mkdir ref_cats
    [gkovacs@lsst-dev03 repo]$ mkdir ref_cats/gaia
    [gkovacs@lsst-dev03 repo]$ mkdir ref_cats/pan-starrs
    [gkovacs@lsst-dev03 repo]$ cd ref_cats/pan-starrs
    [gkovacs@lsst-dev03 pan-starrs]$ tar -xvzf /project/mrawls/ap_verify_ci_hits2015/refcats/ps1_HiTS_2015.tar.gz
    [gkovacs@lsst-dev03 gaia]$ tar -xvzf /project/mrawls/ap_verify_ci_hits2015/refcats/gaia_HiTS_2015.tar.gz 

ProcessCcd
----------

Testcommand:
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ processCcd.py repo/ --rerun proccd_2019-07-03 -C config/processCcd_2019-07-03.py --id visit=419802 ccdnum=1

Collecting visitnums for parallel processing command list:
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ python3 print_visit_nums_fits.py ../hits_cpline_dl/instcals/*.fits.fz | gawk '{print $1}' | uniq | sort -n > ingested_visitnums.txt

Print commands:
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ gawk '{print("processCcd.py repo/ -C config/processCcd_2019-07-03.py --rerun proccd_2019-07-03 --id visit="$1" filter=g")}' ingested_visitnums.txt > procCcd_cmd_2019-07-03

Parallel was not installed on lsst-dev, copied from ubuntu. Fedora version contains weird bibtex/citation cmd line options...
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ xargs -d "\n" parallel -j 6 -- < procCcd_cmd_2019-07-03 |& tee -a procCcd_2019-07-03.log

ImageDifference
---------------

Generate a list of pairs of 2015 visit numbers and 2014 template visit numbers, then the commands.
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ python3 ../pair_calexps_with_template.py ../hits_cpline_dl/instcals/c4d_*fits.fz | sort -n > visit15_template14_pairs.txt
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ gawk '{print("imageDifference.py repo/ -C config/imageDifference_2019-07-04.py --rerun proccd_2019-07-03:imgDiff_2019-07-04 --id visit="$1" filter=g --templateId visit="$2)}' visit15_template14_pairs.txt > imgDiff_cmd_2019-07-04
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ START=`date`; xargs -d "\n" parallel -j 8 -- < imgDiff_cmd_2019-07-04 |& tee -a imgDiff_2019-07-04.log; echo $START >> imgDiff_2019-07-04.log; date >> imgDiff_2019-07-04.log

ApAssociation
--------------
::

    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ make_ppdb.py -c ppdb.db_url=sqlite:///repo/rerun/imgDiff_2019-07-04/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ python3 ../pair_calexps_with_template.py ../hits_cpline_dl/instcals/c4d_*fits.fz | sort -n > visit15_template14_pairs.txt
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ gawk '{ print("python3 run_assoc_2019-07-05.py repo/rerun/imgDiff_2019-07-04 "$1) }' visit15_template14_pairs.txt > assoc_cmd_2019-07-05
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ START=`date`; xargs -d "\n" parallel -j 6 -- < assoc_cmd_2019-07-05 |& tee -a assoc_2019-07-05.log; echo $START >> assoc_2019-07-05.log; date >> assoc_2019-07-05.log

Runing AP association without parallelization: 

::

    [gkovacs@lsst-dev03 imgDiff_2019-07-04]$ mv association.db association_parallel.db
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ make_ppdb.py -c ppdb.db_url=sqlite:///repo/rerun/imgDiff_2019-07-04/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"
    [gkovacs@lsst-dev03 hits_proc2019-07-03]$ START=`date`; source assoc_cmd_2019-07-05 |& tee -a assoc_serial_2019-07-11.log ; echo $START >> assoc_serial_2019-07-11.log; date >> assoc_serial_2019-07-11.log

