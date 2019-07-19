DM-19526 rerun with all Hits2015 visits
----------------------------------------

Create ppdb, generate the ap_pipe commands to process the HiTS15 visits and use their
corresponding HiTS14 calexps as templates. 

2019-07-15
-----------

Now parallelize the 3 fields, they do not overlap so no concurrency issue should occurr.

Pair each 2015 raw file with its 2014 template visit:
::

    [gkovacs@lsst-dev03 hits15_14]$ find ./ -path '*2015*g*fits.fz' | sort | xargs python3  collect_visits_with_template.py > visit15_template14_collection.txt
    
Create a new repo:
::

    [gkovacs@lsst-dev03 hits15_14]$ mkdir ingested/rerun/apPipe_2019-07-15
    [gkovacs@lsst-dev03 hits15_14]$ make_ppdb.py -c ppdb.db_url=sqlite:///ingested/rerun/apPipe_2019-07-15/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"

For whatever reason, the ``--id`` data ids are expanded here to include the ``object`` and ``date`` fields;
if ``--templateId`` specifies a visit id only, the other fields will be inherited from the data specification 
and the butler fails to find the template.
    
::

    [gkovacs@lsst-dev03 hits15_14]$ gawk ' { print("ap_pipe.py ingested/ --calib calibingested/ --rerun apPipe_2019-07-15 -C config/apPipe_2019-06-27.py -c ppdb.db_url=sqlite:///ingested/rerun/apPipe_2019-07-15/association.db -c ppdb.isolation_level=""READ_UNCOMMITTED"" -c ppdb.connection_timeout=240 --id filter=g visit="$1" --templateId visit="$2" object=\""$3"\" date=\"2014-03-03\" & ") }' visit15_template14_collection.txt > ap_pipe_cmd_2019-07-15

    [gkovacs@lsst-dev03 hits15_14]$ START=`date`; source ap_pipe_cmd_2019-07-15 |& tee -a ap_pipe_2019-07-15.log; echo $START >> ap_pipe_2019-07-15.log; date >> ap_pipe_2019-07-15.log


Deleted - 2019-07-13 run
-------------------------

Create ppdb, generate the ap_pipe commands to process one HiTS15 visit and use its
corresponding HiTS14 calexp as template. Note that we need the ``date`` in the templateId
otherwise it inherits the date from ``id`` which won't match.

::

    [gkovacs@lsst-dev03 hits15_14]$ make_ppdb.py -c ppdb.db_url=sqlite:///ingested/rerun/apPipe_2019-07-13/association.db -c ppdb.isolation_level="READ_UNCOMMITTED"
    
    [gkovacs@lsst-dev03 hits15_14]$ gawk ' { print("ap_pipe.py ingested/ --calib calibingested/ --rerun apPipe_2019-07-13 -C config/apPipe_2019-06-27.py -c ppdb.db_url=sqlite:///ingested/rerun/apPipe_2019-07-13/association.db -c ppdb.isolation_level=""READ_UNCOMMITTED"" -c ppdb.connection_timeout=240 --id filter=g visit="$1" --templateId visit="$2" object=\""$3"\" date=\"2014-03-03\" ") }' visit15_template14_pairs.txt > ap_pipe_cmd_2019-07-13

    [gkovacs@lsst-dev03 hits15_14]$ START=`date`; xargs -d "\n" parallel -j 6 -- < ap_pipe_cmd_2019-07-13 |& tee -a ap_pipe_2019-07-13.log; echo $START >> ap_pipe_2019-07-13.log; date >> ap_pipe_2019-07-13.log

    ap_pipe.py ingested/ --calib calibingested/ --rerun apPipe_2019-07-13 -C config/apPipe_2019-06-27.py -c ppdb.db_url=sqlite:///ingested/rerun/apPipe_2019-07-13/association.db -c ppdb.isolation_level=READ_UNCOMMITTED -c ppdb.connection_timeout=240 --id filter=g visit=1 --templateId visit=2 object="3" date="2014-03-03"

No parallel running because of 2014 image proccd may happen concurrently.

::

    [gkovacs@lsst-dev03 hits15_14]$ START=`date`; source ap_pipe_cmd_2019-07-13 |& tee -a ap_pipe_2019-07-13.log; echo $START >> ap_pipe_2019-07-13.log; date >> ap_pipe_2019-07-13.log

This run did not finish in 24 hours, stopped and rerun removed.

