# Running ap_pipe on HSC RC2 data from PDR1

After the new visits from HSC Public Data Release 2 (PDR2) are copied to lsst-dev, we'll recreate an AP-style repo using Public Data Release 1 (PDR1) templates to difference new visits in PDR2. Until then, this repo will demonstrate running `ap_pipe` on a dataset more similar to DRP mode, where the templates are the best seeing images within a data release and all visits in that data release are differenced.

Two repos were generated on lsst-dev:

* `/datasets/hsc/repo/rerun/private/yusra/RC2/ap_pipe_PDR1_bestcase` (HSC-I band)
* `/datasets/hsc/repo/rerun/private/yusra/RC2/ap_pipe_PDR1_bestcase_noSkyCorr` (HSC-G, HSC-R, HSC-I, HSC-Z)

These instructions go into detail on how the second was generated.

To remember the name, the dataset only contains "PDR1" and is "bestcase" in that the templates are constructed from the best seeing images.  The only configuration difference between this rerun and the first was that the HSC-specific skyCorrection was turned off.

Stack version:

```
source /software/lsstsw/stack/loadLSST.bash
setup lsst_distrib -t w_2019_24
```

## Make Templates:


COSMOS visit lists for the PDR1 "RC2" dataset are:

```
export COSMOS_G=11690^11692^11694^11696^11698^11700^11702^11704^11706^11708^11710^11712^29324^29326^29336^29340^29350
export COSMOS_R=1202^1204^1206^1208^1210^1212^1214^1216^1218^1220^23692^23694^23704^23706^23716^23718
export COSMOS_I=1228^1230^1232^1238^1240^1242^1244^1246^1248^19658^19660^19662^19680^19682^19684^19694^19696^19698^19708^19710^19712^30482^30484^30486^30488^30490^30492^30494^30496^30498^30500^30502^30504
export COSMOS_Z=1166^1168^1170^1172^1174^1176^1178^1180^1182^1184^1186^1188^1190^1192^1194^17900^17902^17904^17906^17908^17926^17928^17930^17932^17934^17944^17946^17948^17950^17952^17962
export COSMOS_Y=318^322^324^326^328^330^332^344^346^348^350^352^354^356^358^360^362^1868^1870^1872^1874^1876^1880^1882^11718^11720^11722^11724^11726^11728^11730^11732^11734^11736^11738^11740^22602^22604^22606^n^22626^22628^22630^22632^22642^22644^22646^22648^22658^22660^22662^22664
```

Manually choose visits with the top third best seeing.
```
export COSMOS_G_SEEING=11704^11702^11694^11698^11692^11696
export COSMOS_R_SEEING=23694^23692^1208^1218^23704
export COSMOS_I_SEEING=19682^19680^19684^19662^19658^30502^30490^30500^19696^19694
export COSMOS_Z_SEEING=17928^17908^17906^17926^17904^17902^17930^17962^17946^17944
export COSMOS_Y_SEEING=22602^22604^22644^22632^22628^22642^22626^22630^22664^11730^22662^22606^11736^22646^22660^11718^11732
```

These visit lists are sorted by ascending visit FHWM, so reducing the number of visits can be accomplished by removing visits from the end.


Generate templates without skyCorrection for these best seeing visits:

```
coaddDriver.py  /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_templates_noSkyCorr --batch-type=slurm --mpiexec='-bind-to socket' --job coaddCG --time 600 --cores 24  --id tract=9813 filter=HSC-G --selectId ccd=0..8^10..103 visit=$COSMOS_G_SEEING  --config makeCoaddTempExp.doApplySkyCorr=False --no-versions
coaddDriver.py  /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_templates_noSkyCorr --batch-type=slurm --mpiexec='-bind-to socket' --job coaddCR --time 600 --cores 24  --id tract=9813 filter=HSC-R --selectId ccd=0..8^10..103 visit=$COSMOS_R_SEEING  --config makeCoaddTempExp.doApplySkyCorr=False --no-versions
coaddDriver.py  /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_templates_noSkyCorr --batch-type=slurm --mpiexec='-bind-to socket' --job coaddCI --time 600 --cores 24  --id tract=9813 filter=HSC-I --selectId ccd=0..8^10..103 visit=$COSMOS_I_SEEING  --config makeCoaddTempExp.doApplySkyCorr=False --no-versions
coaddDriver.py  /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_templates_noSkyCorr --batch-type=slurm --mpiexec='-bind-to socket' --job coaddCZ --time 600 --cores 24  --id tract=9813 filter=HSC-Z --selectId ccd=0..8^10..103 visit=$COSMOS_Z_SEEING  --config makeCoaddTempExp.doApplySkyCorr=False --no-versions
coaddDriver.py  /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_templates_noSkyCorr --batch-type=slurm --mpiexec='-bind-to socket' --job coaddCY --time 600 --cores 24  --id tract=9813 filter=HSC-Y --selectId ccd=0..8^10..103 visit=$COSMOS_Y_SEEING  --config makeCoaddTempExp.doApplySkyCorr=False --no-versions
```


## Set up the ppdb:

```
mkdir /datasets/hsc/repo/rerun/private/yusra/RC2/ap_pipe_PDR1_bestcase_noSkyCorr
make_ppdb.py --configfile configfile.py
```

where `configfile.py` contains:
```
$ cat configfile.py
config.ppdb.db_url="sqlite:////datasets/hsc/repo/rerun/private/yusra/RC2/ap_pipe_PDR1_bestcase_noSkyCorr/association.db"
config.ppdb.isolation_level="READ_UNCOMMITTED"
config.ppdb.connection_timeout=540
config.differencer.doWriteMatchedExp=True
```

Note that only the first two lines are used by `make_ppdb.py` the 3rd and 4th line are used by`ap_pipe`:

## Run ap_pipe:

Out of habit, I use  Hsin-Fang's latest monthly rerun `RC/w_2019_22/DM-19244` as input.  This will just redirect up the parent tree to the raws in `/datasets/hsc` when we ask for raws.

```
ap_pipe.py /datasets/hsc/repo --calib /datasets/hsc/repo/CALIB/ --template /datasets/hsc/repo/rerun/private/yusra/RC2/ap_pipe_templates_noSkyCorr --rerun RC/w_2019_22/DM-19244:private/yusra/RC2/ap_pipe_PDR1_bestcase_noSkyCorr --configfile ./configfile.py --no-versions --longlog --id visit=$COSMOS_I ccd=0..8^10..103
```

Note that `ccd=0..8^10..103`. We generally skip ccd=9, since more than half the chip has no usable data.

In practice launched these as slurm jobs for all visits in `$COSMOS_G`, `$COSMOS_R`, `$COSMOS_I`, and `$COSMOS_Z`.



# Notes:

* If I launched more than approx 30 processes running at the same time I got sqlite connection errors. This felt quite limiting.
```
sqlite3.OperationalError: database is locked
sqlalchemy.exc.OperationalError: (sqlite3.OperationalError) database is locked
```

* As the database grew in size, `AssociationTask` started taking 15-20 minutes to run per ccd.
* `AssociationTask` does not log.
