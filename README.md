# ap_pipe-notebooks
A repository for notebooks and other useful code intended to  analyze 
and troublueshoot various ap_pipe processing runs.

**data**: Known variables in the HiTS 2015 dataset, from the HiTS DR1

**false_positives**: Work during the June 2019 sprint to characterize false positives in difference imaging

*Note:* [RFC-642](https://jira.lsstcorp.org/browse/RFC-642) renamed PPDB
into APDB, many module names and method names have changed as a result.
Python modules in this package have been updated for that change but notebooks
have not been changed. This can potentially break some or all notebooks,
be prepared to update them.

*Note:* [DM-34627](https://jira.lsstcorp.org/browse/DM-34627) moved the functions in apdbPlots.py,
coaddAnalysis.py, diaObjectAnalysis.py, and plotLightCurve.py into
legacyApdbUtils.py, legacyCoaddAnalysis.py, and legacyPlotUtils.py in analysis_ap.
In order to get older notebooks working with these new functions instead of
the older ones, change the import calls to point to the new function in analysis_ap,
e.g. "from lsst.analysis.ap import legacyApbUtils, legacyCoaddAnalysis, legacyPlotUtils".
