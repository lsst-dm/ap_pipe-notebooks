import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import functools
import operator
from astropy import units as u

import lsst.daf.persistence as dafPersist
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay
import lsst.geom
from lsst.ap.association import UnpackApdbFlags, TransformDiaSourceCatalogConfig
import lsst.afw.cameraGeom as cameraGeom
from lsst.obs.decam import DarkEnergyCamera

from diaObjectAnalysis import loadAllApdbObjects, loadAllApdbSources
from plotLightcurve import getTemplateCutoutGen2 as getTemplateCutout
"""
Collection of plots that can be made using info in the APDB.

In some cases, additional inputs may be needed, and that will be
clearly indicated in the docstrings.

Future plots could include:
- plot DIA Sources on focal plane on one visit or ccd at a time
- plot various totFlux vs psFlux and psFlux vs apFlux (or magnitudes)
- plot stats about the images comprising the template
"""


def addVisitCcdToSrcTable(sourceTable, instrument='DECam', gen='gen2', butler=None):
    """Add visit and ccd columns to sourceTable dataframe.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    instrument : `str`, optional
        Defaults to 'DECam'; also supports gen2 'HSC' or any gen3 instrument.
    gen : `str`, optional
        Either 'gen2' or 'gen3'.
    butler : `lsst.daf.butler.Butler` or None (gen2), optional
        Butler in the repository corresponding to the output of an ap_pipe run.

    Returns
    -------
    sourceTable : `pandas.core.frame.DataFrame`
        The same as the input sourceTable, with new visit and ccd columns.
    """
    if instrument == 'DECam' and gen == 'gen2':
        sourceTable['ccd'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[-2:])
        sourceTable['visit'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[:-2])
    elif instrument == 'HSC' and gen == 'gen2':
        print('WARNING: your visit values may be suspect')
        sourceTable['visit'] = sourceTable.ccdVisitId.apply(lambda x: int(x/200.))
        sourceTable['ccd'] = sourceTable.ccdVisitId.apply(lambda x: int(np.round((x/200. - int(x/200.))*200)))
    elif gen == 'gen3':  # should work for all instruments
        instrumentDataId = butler.registry.expandDataId(instrument=instrument)
        packer = butler.registry.dimensions.makePacker("visit_detector", instrumentDataId)
        dataId = packer.unpack(sourceTable.ccdVisitId)
        sourceTable['visit'] = dataId['visit']
        sourceTable['ccd'] = dataId['detector']
    elif gen == 'gen2' and instrument not in ['DECam', 'HSC']:
        raise ValueError('Keyword `instrument` is case sensitive. "DECam" and "HSC" work with gen2.')
    return sourceTable


def plot2axes(x1, y1, x2, y2, xlabel, ylabel1, ylabel2, y1err=None, y2err=None, title=''):
    """Generic plot framework for showing one y-variable on the left
    and another on the right. The x-axis is shared on the bottom.

    Parameters
    ----------
    x1 : `list` or `array`
    y1 : `list` or `array`
    x2 : `list` or `array`
    y2 : `list` or `array`
    xlabel : `str`
        Label for shared x axis
    ylabel1 : `str`
        Label for y1 axis
    ylabel2 : `str`
        Label for y2 axis
    title : `str`
        Title for the plot, optional.
    """
    fig, ax1 = plt.subplots(figsize=(9, 3))
    color = 'C1'
    if y1err is not None:
        ax1.errorbar(x1, y1, yerr=y1err, color=color, marker='o', ls=':')
    else:
        ax1.plot(x1, y1, color=color, marker='o', ls=':')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    for label in ax1.get_xticklabels():
        label.set_ha("right")
        label.set_rotation(45)
        label.set_size("smaller")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'C0'
    ax2.set_ylabel(ylabel2, color=color)  # we already handled the x-label above
    if y2err is not None:
        ax2.errorbar(x2, y2, yerr=y2err, color=color, marker='o', ls=':')
    else:
        ax2.plot(x2, y2, color=color, marker='o', ls=':')
    ax2.tick_params(axis='y', labelcolor=color)
    plt.title(title)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped


def plotDiaSourcesPerVisit(repo, sourceTable, title='', gen='gen2', instrument='DECam', collections=[]):
    """Plot DIA Sources per visit.

    The plot will have two y-axes: number of DIA Sources per square degree and
    median FWHM per ixx or iyy in pixels.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    title : `str`
        Title for the plot, optional.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly
    """
    ccdArea, visitArea = getCcdAndVisitSizeOnSky(repo, sourceTable, gen, instrument, collections)
    traceRadius = np.sqrt(0.5) * np.sqrt(sourceTable.ixxPSF + sourceTable.iyyPSF)
    sourceTable['seeing'] = 2*np.sqrt(2*np.log(2)) * traceRadius
    visitGroup = sourceTable.groupby('visit')
    plot2axes(visitGroup.visit.first().values,
              visitGroup.ccd.count().values/visitArea,
              visitGroup.visit.first().values,
              visitGroup.seeing.median().values,
              'Visit',
              'Number of DIA Sources (per sq. deg.)',
              'Median FWHM per ixx/iyy (pixels)',
              title=title)


def plotDiaSourcesPerNight(sourceTable, title=''):
    """Plot DIA Sources per night.

    The plot will have two y-axes: mean number of DIA Sources per visit and
    number of visits per night.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
        NOT a view into or slice of a dataframe!
    title : `str`
        Title for the plot, optional.
    """
    date_times = pd.to_datetime(sourceTable['midPointTai'],
                                unit='D',
                                origin=pd.Timestamp('1858-11-17'))
    sourceTable['date'] = date_times.dt.date
    night_count = sourceTable.groupby(['date', 'visit']).count()
    visits_per_night = night_count.groupby('date').count()
    pervisit_per_night = night_count.groupby('date').mean()
    pervisit_per_night_std = night_count.groupby('date').std()
    pervisit_per_night_err = pervisit_per_night_std['x']/np.sqrt(visits_per_night['x'])
    plot2axes(visits_per_night.index,
              pervisit_per_night['x'],
              visits_per_night.index,
              visits_per_night['x'],
              'Night',
              'Mean DIA Sources per visit',
              'Number of Visits',
              y1err=pervisit_per_night_err,
              title=title)
    plt.ylim(0, np.max(visits_per_night['x'].values + 1))


def ccd2focalPlane(x, y, ccd, camera):
    """Retrieve focal plane coordinates.

    Parameters
    ----------
    x : `int` or `float`
        X-coordinate from ccd bbox.
    y : `int` or `float`
        Y-coordinate from ccd bbox.
    ccd : `int`, or can be cast as int
        The ccd being considered.
    camera : `lsst.afw.cameraGeom.Camera`
        Camera obtained from the butler.

    Returns
    -------
    point[0] : `int` or `float`
        X-coordinate in focal plane.
    point[1] : `int` or `float`
        Y-coordinate in focal plane.
    """
    detector = camera[int(ccd)]
    point = detector.transform(lsst.geom.Point2D(x, y),
                               cameraGeom.PIXELS, cameraGeom.FOCAL_PLANE)
    return point[0], point[1]


def getCcdCorners(butler, sourceTable, gen='gen2', instrument='DECam', collections=[]):
    """Get corner coordinates for a range of ccds.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler` or `lsst.daf.butler.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly

    Returns
    -------
    corners : `pandas.core.frame.DataFrame`
        Dataframe containing focal plane coordinates for all the ccd corners.
    """
    cornerList = []
    visits = np.unique(sourceTable['visit'])
    ccds = np.unique(sourceTable['ccd'])
    ccdMin = int(np.min(ccds))
    ccdMax = int(np.max(ccds))
    visit = int(visits[0])  # shouldn't matter what visit we use
    for ccd in range(ccdMin, ccdMax+1):
        try:
            if gen == 'gen2':
                bbox = butler.get('calexp_bbox', dataId={visit: visit, ccd: ccd})
            else:  # gen3
                bbox = butler.get('calexp.bbox', collections=collections,
                                  instrument=instrument, visit=visit, detector=ccd)
        except (dafPersist.NoResults, LookupError):  # silently skip any ccds that don't exist
            for visit in visits[1:]:  # try the other visits just in case
                visit = int(visit)
                try:
                    bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
                except (dafPersist.NoResults, LookupError):
                    continue
                break
        else:
            if gen == 'gen2':
                camera = butler.get('camera')
            else:  # gen3
                if instrument == 'DECam':
                    camera = DarkEnergyCamera().getCamera()
                else:
                    raise NotImplementedError
            cornerList.append([visit, ccd] + [val for pair in bbox.getCorners()
                              for val in ccd2focalPlane(pair[0], pair[1], ccd, camera)])
    corners = pd.DataFrame(cornerList)
    corners['width'] = corners[6] - corners[8]
    corners['height'] = corners[7] - corners[5]
    return corners


def getCcdAndVisitSizeOnSky(repo, sourceTable, gen='gen2', instrument='DECam', collections=[],
                            visit=None, detector=None):
    """Estimate the area of one CCD and one visit on the sky, in square degrees.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly
    visit : `int` or None, optional
        Specific visit to use when loading representative calexp.
    detector : `int` or None, optional
        Specific detector (ccd) to use when loading representative calexp.

    Returns
    -------
    ccdArea : `float`
        Area covered by one detector (CCD) on the sky, in square degrees
    visitArea :
        Area covered by a visit with all detectors (CCDs) on the sky, in square degrees
    """
    visits = np.unique(sourceTable.visit)
    ccds = np.unique(sourceTable.ccd)
    nGoodCcds = len(ccds)
    if gen == 'gen2':
        butler = dafPersist.Butler(repo)
        if visit is None:
            visit = int(visits[0])
        if detector is None:
            ccd = int(ccds[0])
        calexp = butler.get('calexp', visit=visit, ccd=ccd)
        bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
    else:  # gen3
        butler = dafButler.Butler(repo)
        if visit is None:
            visit = int(visits[0])
        if detector is None:
            detector = int(ccds[0])
        calexp = butler.get('calexp', collections=collections,
                            instrument=instrument, visit=visit, detector=detector)
        bbox = butler.get('calexp.bbox', collections=collections,
                          instrument=instrument, visit=visit, detector=detector)
    pixelScale = calexp.getWcs().getPixelScale().asArcseconds()
    ccdArea = (pixelScale*pixelScale*bbox.getArea()*u.arcsec**2).to(u.deg**2).value
    visitArea = ccdArea * nGoodCcds
    return ccdArea, visitArea


def plotDiaSourceDensityInFocalPlane(repo, sourceTable, cmap=mpl.cm.Blues, title='',
                                     gen='gen2', instrument='DECam', collections=[]):
    """Plot average density of DIA Sources in the focal plane (per CCD).

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    cmap : `matplotlib.colors.ListedColormap`
        Matplotlib colormap.
    title : `str`
        String to append to the plot title, optional.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly
    """
    ccdArea, visitArea = getCcdAndVisitSizeOnSky(repo, sourceTable, gen, instrument, collections)
    nVisits = len(np.unique(sourceTable['visit'].values))
    ccdGroup = sourceTable.groupby('ccd')
    ccdSourceCount = ccdGroup.visit.count().values/nVisits/ccdArea
    # DIA Source count per visit per square degree, for each CCD
    if gen == 'gen2':
        butler = dafPersist.Butler(repo)
    else:  # gen3
        butler = dafButler.Butler(repo)
    corners = getCcdCorners(butler, sourceTable, gen, instrument, collections)
    norm = mpl.colors.Normalize(vmin=np.min(ccdSourceCount), vmax=np.max(ccdSourceCount))

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111, aspect='equal')

    for index, row in corners.iterrows():
        try:
            averageFocalPlane = ccdGroup.get_group(int(row[1])).x.count()/nVisits/ccdArea
        except KeyError:
            averageFocalPlane = 0  # plot normalization will be weird but it won't fall over
        ax1.add_patch(patches.Rectangle((row[7], row[6]),
                      -row.height,
                      -row.width,
                      fill=True,
                      color=cmap(norm(averageFocalPlane))))
        ax1.text(row[7]-row.height/2, row[6]-row.width/2, '%d' % (row[1]), fontsize=12)
        plt.plot(row[7]-row.height/2, row[6]-row.width/2, ',')
    ax1.set_title('Mean DIA Source density in focal plane coordinates %s' % (title))
    ax1.set_xlabel('Focal Plane X', size=16)
    ax1.set_ylabel('Focal Plane Y', size=16)
    ax1 = plt.gca()
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, fraction=0.04, pad=0.04)
    cb.set_label('DIA Sources per sq. deg.', rotation=90)


def loadTables(repo, dbName='association.db', isVerify=False, dbType='sqlite',
               badFlagList=['base_PixelFlags_flag_bad',
                            'base_PixelFlags_flag_suspect',
                            'base_PixelFlags_flag_saturatedCenter'],
               instrument='DECam', gen='gen2', schema=None):
    """Load DIA Object and DIA Source tables from an APDB.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    dbName : `str`
        Name of APDB.
    isVerify : `bool`
        Is this an ap_verify run instead of an ap_pipe run?
        If True, the APDB is one level above repo on disk.
        If False, the APDB is in repo (default).
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    badFlagList :  `list`
        Names of flags presumed to each indicate a DIA Source is garbage.
    instrument : `str`, one of either 'DECam' or 'HSC', default is 'DECam'
        Needed to properly add the "ccd" and "visit" columns to the sourceTable,
        and for all things gen3
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    schema : `str`, optional
        Required if dbType is postgres

    Returns
    -------
    objTable : `pandas.core.frame.DataFrame`
        DIA Object table loaded from the APDB.
    srcTable : `pandas.core.frame.DataFrame`
        DIA Source table loaded from the APDB.
    goodObj : `pandas.core.frame.DataFrame`
        A subset of objTable containing only DIA Objects composed entirely of good DIA Sources.
    goodSrc : `pandas.core.frame.DataFrame`
        A subset of srcTable containing only good DIA Sources.
    """
    if dbType == 'sqlite' and not isVerify:  # APDB is in repo (default)
        dbPath = os.path.abspath(os.path.join(repo, dbName))
    elif dbType == 'sqlite' and isVerify:  # APDB is one level above repo (ap_verify sqlite case)
        repoUpOne = os.path.dirname(repo)
        dbPath = os.path.abspath(os.path.join(repoUpOne, dbName))
    elif dbType == 'postgres':
        dbPath = dbName
    else:
        raise ValueError('database type not understood')

    if gen == 'gen3':
        butler = dafButler.Butler(repo)
    else:
        butler = None

    objTable = loadAllApdbObjects(dbPath, dbType=dbType, schema=schema)
    srcTable = loadAllApdbSources(dbPath, dbType=dbType, schema=schema)
    srcTable = addVisitCcdToSrcTable(srcTable, instrument=instrument, gen=gen, butler=butler)
    flagTable, srcTableFlags, flagFilter, \
        goodSrc, goodObj = makeSrcTableFlags(srcTable, objTable, badFlagList=badFlagList,
                                             gen=gen, instrument=instrument, repo=repo)
    return objTable, srcTable, goodObj, goodSrc


def makeSrcTableFlags(sourceTable, objectTable,
                      badFlagList=['base_PixelFlags_flag_bad',
                                   'base_PixelFlags_flag_suspect',
                                   'base_PixelFlags_flag_saturatedCenter'],
                      gen='gen2', instrument='DECam', repo=None):
    """Apply flag filters to a DIA Source and a DIA Object table.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    objectTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Objects from the same APDB.
    badFlagList :  `list`
        Names of flags presumed to each indicate a DIA Source is garbage.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    repo : `str`, optional
        Repository in which to load a butler, used with gen3 only

    Returns
    -------
    flagTable : `pandas.core.frame.DataFrame`
        Dataframe containing unpacked DIA Source flag values.
    sourceTableFlags : `pandas.core.frame.DataFrame`
        Dataframe resulting from from merging sourceTable with flagTable.
    flagFilter : `pandas.core.series.Series` of `bool`
        Single column of booleans of length len(sourceTable).
        The value of flagFilter is True if one or more flags in badFlagList is True.
    goodSrc : `pandas.core.frame.DataFrame`
        Dataframe containing only DIA Sources from sourceTable with no bad flags.
    goodObj : `pandas.core.frame.DataFrame`
        Dataframe containing only DIA Objects from objectTable entirely composed
        of DIA Sources with no bad flags.
    """
    if gen == 'gen3':
        butler = dafButler.Butler(repo)
    else:
        butler = None
    sourceTable = addVisitCcdToSrcTable(sourceTable, instrument=instrument, gen=gen, butler=butler)
    config = TransformDiaSourceCatalogConfig()
    unpacker = UnpackApdbFlags(config.flagMap, 'DiaSource')
    flagValues = unpacker.unpack(sourceTable['flags'], 'flags')
    flagTable = pd.DataFrame(flagValues, index=sourceTable.index)
    sourceTableFlags = pd.merge(sourceTable, flagTable, left_index=True, right_index=True)
    badFlags = [sourceTableFlags[flag] for flag in badFlagList]
    flagFilter = functools.reduce(operator.or_, badFlags)
    noFlagFilter = ~flagFilter
    goodSrc = sourceTableFlags.loc[noFlagFilter]
    goodObjIds = set(sourceTableFlags.loc[noFlagFilter, 'diaObjectId'])
    goodObj = objectTable.loc[objectTable['diaObjectId'].isin(goodObjIds)]
    return flagTable, sourceTableFlags, flagFilter, goodSrc, goodObj


def plotHitsSourcesOnSky(sourceTable, title=''):
    """Plot DIA Sources from three DECam HiTS fields on the sky.

    Can also be used to plot DIA Objects instead of DIA Sources.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    title : `str`
        Title for the plot, optional.
    """
    plt.figure(figsize=(9, 7))
    ax1 = plt.subplot2grid((100, 100), (0, 55),
                           rowspan=50, colspan=45)  # 1 single HiTS field, on the right
    ax2 = plt.subplot2grid((100, 100), (0, 0),
                           rowspan=90, colspan=50)  # 2 overlapping HiTS fields, on the left

    ax1Filter = (sourceTable['decl'] > -2)
    ax2Filter = (~ax1Filter)

    ax1.scatter(sourceTable.loc[ax1Filter, 'ra'],
                sourceTable.loc[ax1Filter, 'decl'],
                c='C0', s=2, alpha=0.2)
    ax2.scatter(sourceTable.loc[ax2Filter, 'ra'],
                sourceTable.loc[ax2Filter, 'decl'],
                c='C0', s=2, alpha=0.2)

    ax1.set_xlabel('RA (deg)')
    ax2.set_xlabel('RA (deg)')
    ax1.set_ylabel('Dec (deg)')
    ax2.set_ylabel('Dec (deg)')

    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.invert_xaxis()
    ax2.invert_xaxis()

    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position('right')

    plt.title(title)
    plt.subplots_adjust(wspace=0.1, hspace=0)


def plotSourceLocationsOnObject(objId, objTable, sourceTable, repo, goodSrc=None, plotIdx=1):
    """Plot a difference imaging template used for a DIA Object in a postage
    stamp cutout, along with the locations of the DIA Object and all
    constituent DIA Sources.

    This is a deprecated gen2 function for the time being.

    Parameters
    ----------
    objId : `int`
        DIA Object ID.
    objTable : `pandas.core.frame.DataFrame`
        DIA Object table loaded from the APDB.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB, plus unpacked flag values.
        This is an output called "sourceTableFlags" from the makeSrcTableFlags function.
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    goodSrc : `pandas.core.frame.DataFrame` or `None`
        Filtered subset of sourceTable. Sources in this table will be plotted
        in green, while those not in this table will be plotted in red.
        If `None`, all sources will be plotted in blue.
    plotIdx : `int`
        Frame number to put the resulting plot in, default 1.

    Notes
    -----
    The sourceTable is assumed to have flag fields added to it already!
    """
    butler = dafPersist.Butler(repo)
    afwDisplay.setDefaultBackend('matplotlib')

    objRa = objTable.loc[objTable['diaObjectId'] == objId, 'ra'].values[0]
    objDec = objTable.loc[objTable['diaObjectId'] == objId, 'decl'].values[0]
    srcIds = list(sourceTable.loc[sourceTable['diaObjectId'] == objId, 'diaSourceId'].values)
    srcRas = sourceTable.loc[sourceTable['diaObjectId'] == objId, 'ra'].values
    srcDecs = sourceTable.loc[sourceTable['diaObjectId'] == objId, 'decl'].values
    srcVisits = sourceTable.loc[sourceTable['diaObjectId'] == objId, 'visit'].values
    srcCcds = sourceTable.loc[sourceTable['diaObjectId'] == objId, 'ccd'].values
    dataIds = [{'visit': int(visit), 'ccd': int(ccd)} for (visit, ccd) in zip(srcVisits, srcCcds)]

    objLoc = lsst.geom.SpherePoint(objRa, objDec, lsst.geom.degrees)
    calexpFirst = butler.get('calexp', dataId=dataIds[0])
    templateCutout = getTemplateCutout(calexpFirst, repo, objLoc)
    if templateCutout is not None:
        imageToShow = templateCutout
    else:
        imageToShow = calexpFirst.getCutout(objLoc, size=lsst.geom.Extent2I(30, 30))
        print('Failed to retrieve template; image displayed is first processed image')
    wcs = imageToShow.getWcs()

    disp = afwDisplay.Display(plotIdx, reopenPlot=True)
    disp.setMaskTransparency(100)
    disp.scale('asinh', 'zscale', Q=4)
    disp.mtv(imageToShow, title=objId)

    with disp.Buffering():  # obviously
        coordObj = wcs.skyToPixel(objLoc)
        disp.dot('+', *coordObj, ctype='C0', size=5)
        for dataId, srcId, ra, dec in zip(dataIds, srcIds, srcRas, srcDecs):
            calexp = butler.get('calexp', dataId=dataId)
            psf = calexp.getPsf()
            psfSize = psf.computeShape().getDeterminantRadius()*2.355  # sigma to FWHM
            coordSrc = wcs.skyToPixel(lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees))
            if goodSrc is not None:
                if srcId in goodSrc['diaSourceId'].values:
                    disp.dot('o', *coordSrc, ctype='C3', size=psfSize)  # green
                else:
                    disp.dot('o', *coordSrc, ctype='C2', size=psfSize)  # red
            else:
                disp.dot('o', *coordSrc, ctype='C0', size=psfSize)  # blue


def plotSeeingHistogram(repo, sourceTable, ccd=35, gen='gen2', instrument='DECam', collections=[]):
    """Plot distribution of visit seeing.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    ccd : `int`
        The ccd being considered, default 35.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly
    """
    fwhm = pd.DataFrame()
    visits = np.unique(sourceTable['visit'])
    radii = []
    if gen == 'gen2':
        butler = dafPersist.Butler(repo)
        for visit in visits:
            calexp = butler.get('calexp_sub', visit=int(visit),
                                ccd=ccd, bbox=lsst.geom.Box2I())
            psf = calexp.getPsf()
            psfSize = psf.computeShape().getDeterminantRadius()
            radii.append(psfSize*2.355)  # convert sigma to FWHM
    else:  # gen3
        butler = dafButler.Butler(repo)
        for visit in visits:
            psf = butler.get('calexp.psf', instrument=instrument,
                             collections=collections,
                             visit=int(visit), detector=ccd)
            psfSize = psf.computeShape().getDeterminantRadius()
            radii.append(psfSize*2.355)  # convert sigma to FWHM
        # Get just one calexp for WCS purposes
        calexp = butler.get('calexp', collections=collections,
                            instrument=instrument,
                            visit=int(visit), detector=ccd)
    fwhm['visit'] = pd.Series(visits)
    fwhm['radius'] = pd.Series(radii, index=fwhm.index)
    pixelScale = calexp.getWcs().getPixelScale().asArcseconds()  # same for all visits
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.hist(fwhm['radius'].values, alpha=0.5)
    plt.xlabel('Seeing FWHM (pixels)')
    plt.ylabel('Visit count')
    secax = ax.secondary_xaxis('top', functions=(lambda x: x*pixelScale, lambda x: x/pixelScale))
    secax.set_xlabel('Seeing FWHM (arcseconds)')


def plotDiaSourcesInFocalPlane(repo, sourceTable, gridsize=(400, 400), title='',
                               gen='gen2', instrument='DECam', collections=[]):
    """Plot DIA Source locations in the focal plane.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    gridsize : `tuple` of form (int, int)
        Number of hexagons in the (x, y) directions for the hexbin plot.
    title : `str`
        String to append to the plot title, optional.
    gen : `str`, optional
        Either 'gen2' or 'gen3'
    instrument : `str`, optional
        Default is 'DECam', used with gen3 butler only
    collections : `list` or `str`, optional
        Must be provided for gen3 to load the camera properly
    """
    if gen == 'gen2':
        butler = dafPersist.Butler(repo)
        camera = butler.get('camera')
    else:
        butler = dafButler.Butler(repo)
        if instrument == 'DECam':
            camera = DarkEnergyCamera().getCamera()
        else:
            raise NotImplementedError
    corners = getCcdCorners(butler, sourceTable, gen, instrument, collections)
    xFP_list = []
    yFP_list = []
    for index, row in sourceTable.iterrows():
        xFP, yFP = ccd2focalPlane(row['x'], row['y'], row['ccd'], camera=camera)
        xFP_list.append(xFP)
        yFP_list.append(yFP)
    xFP_Series = pd.Series(xFP_list, index=sourceTable.index)
    yFP_Series = pd.Series(yFP_list, index=sourceTable.index)

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111, aspect='equal')

    for index, row in corners.iterrows():
        ax1.add_patch(patches.Rectangle((row[7], row[6]),
                                        -row.height,
                                        -row.width,
                                        fill=False))
        ax1.text(row[7] - row.height/2, row[6] - row.width/2, '%d' % (row[1]))
        plt.plot(row[7] - row.height/2, row[6] - row.width/2, ',')

    # somehow x and y are switched... geometry is hard
    ax1.hexbin(yFP_Series, xFP_Series, gridsize=gridsize, bins='log', cmap='Blues')
    ax1.set_title('DIA Sources in focal plane coordinates %s' % (title))

    ax1.set_xlabel('Focal Plane X', size=16)
    ax1.set_ylabel('Focal Plane Y', size=16)
    ax1.invert_yaxis()
    ax1.invert_xaxis()


def plotDiaSourcesOnSkyGrid(repo, sourceTable, title=None, color='C0', size=0.1):
    """Make a multi-panel plot of DIA Sources for each visit on the sky.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    title : `str`
        String for overall figure title, optional.
    color : `str`
        Color to use for the scatter plot, optional (default is C0 blue).
    size : `float`
        Size for points with marker style '.'.
    """
    visits = np.unique(sourceTable['visit'])
    nVisits = len(visits)
    if np.floor(np.sqrt(nVisits)) - np.sqrt(nVisits) == 0:
        squareGridSize = np.int(np.sqrt(nVisits))
    else:
        squareGridSize = np.int(np.sqrt(nVisits)) + 1
    fig = plt.figure(figsize=(9, 9))
    for count, visit in enumerate(np.unique(sourceTable['visit'].values)):
        idx = sourceTable.visit == visit
        ax = fig.add_subplot(squareGridSize, squareGridSize, count + 1, aspect='equal')
        ax.scatter(sourceTable.ra[idx], sourceTable.decl[idx], c=color,
                   marker='.', s=size, alpha=0.2)
        ax.set_title(visit, size=8)
        ax.invert_xaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.subplots_adjust(wspace=0)
    if title:
        fig.suptitle(title)


def plotFlagHist(sourceTable, title=None,
                 badFlagList=['base_PixelFlags_flag_bad',
                              'base_PixelFlags_flag_suspect',
                              'base_PixelFlags_flag_saturatedCenter']):
    """Plot a histogram of how often each pixel flag occurs in DIA Sources.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    title : `str`
        String for overall figure title, optional.
    badFlagList : `list`, optional
        Flag names to plot in red, presumed to indicate a DIA Source is garbage.
    """
    config = TransformDiaSourceCatalogConfig()
    unpacker = UnpackApdbFlags(config.flagMap, 'DiaSource')
    flagValues = unpacker.unpack(sourceTable['flags'], 'flags')
    labels = list(flagValues.dtype.names)
    flagTable = pd.DataFrame(flagValues, index=sourceTable.index)
    flagSum = flagTable.sum()
    flagsToPlot = [count for count in flagSum.values]
    assert len(flagsToPlot) == len(labels)

    flagColors = []
    for label in labels:
        if label in badFlagList:
            flagColors.append('C3')
        else:
            flagColors.append('C0')

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.barh(labels, flagsToPlot, color=flagColors)
    fig.subplots_adjust(left=0.35)
    ax.set_xlabel('Number of flagged DIASources')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title(title)


def plotFluxHistSrc(srcTable1, srcTable2=None, fluxType='psFlux',
                    label1=None, label2=None, title=None, ylog=False,
                    color1='#2979C1', color2='#Bee7F5',
                    binmin=-3e3, binmax=3e3, nbins=200):
    """Plot distribution of fluxes from 1-2 DIASource tables.

    Parameters
    ----------
    srcTable1 : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    srcTable2 : `pandas.core.frame.DataFrame`, optional
        Second pandas dataframe with DIA Sources from an APDB.
    fluxType : `str`, optional
        Choose from psFlux (default), totFlux, or apFlux.
    label1 : `str`, optional
        Label for srcTable1.
    label2 : `str`, optional
        Label for srcTable2.
    title : `str`, optional
        Plot title.
    ylog : `bool`, optional
        Plot the y-axis on a log scale? Default False
    color1 : `str`, optional
        Color for srcTable1.
    color2 : `str`, optional
        Color for srcTable2.
    binmin : `float`, optional
        Minimum x-value for bins.
    binmax : `float`, optional
        Maximum x-value for bins.
    nbins : `int`, optional
        Number of histogram bins.

    """
    plt.figure()
    plt.xlabel(fluxType, size=12)
    plt.ylabel('DIA Source count', size=12)
    bins = np.linspace(binmin, binmax, nbins)
    if ylog:
        plt.yscale('log')
    plt.hist(srcTable1[fluxType].values, bins=bins, color=color1, label=label1)
    if srcTable2 is not None:
        plt.hist(srcTable2[fluxType].values, bins=bins, color=color2, label=label2)
    if label1:
        plt.legend(frameon=False, fontsize=12)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title(title)


def source_magnitude_histogram(repo, sourceTable, bandToPlot, instrument,
                               collections, detectorToUse=42, binsize=0.2,
                               min_magnitude=None, max_magnitude=None,
                               badFlagList=['base_PixelFlags_flag_bad',
                                            'base_PixelFlags_flag_suspect',
                                            'base_PixelFlags_flag_saturatedCenter']):
    """Plot magnitudes of sources from the source catalog for DIA Sources.

    This is a gen3-only function!

    Note that the values plotted are from the 'src' (catalog corresponding to
    the processed visit image), not the APDB or difference image source catalog.

    Parameters
    ----------
    repo : `str`
        Gen3 repository containing 'collections'
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    bandToPlot : `str`
        Typically one of 'g', 'r', 'i', 'z', or 'y'
    instrument : `str`
        e.g., 'HSC' or 'DECam'
    collections : `str`
        Collection within gen3 'repo' to use
    detectorToUse : `int`, optional
        Detector to use for all the plots, default 42
    binsize : `float`, optional
        Histogram bin size, default 0.2
    min_magnitude : `float` or None, optional
        Set plot x-axis minimum
    max_magnitude : `float` or None, optional
        Set plot x-axis maximum
    badFlagList : `list`, optional
        Exclude sources with flags in this list.
    """

    visits = np.unique(sourceTable['visit'])
    area, _ = getCcdAndVisitSizeOnSky(repo, sourceTable, gen='gen3', detector=detectorToUse,
                                      instrument=instrument, collections=collections)
    butler = dafButler.Butler(repo)

    fig, ax = plt.subplots(figsize=(8, 5))
    for visit in visits:
        band = sourceTable.loc[sourceTable['visit'] == visit, 'filterName'].values[0]
        if band != bandToPlot:
            continue
        if band == 'g':
            color = 'C2'
        elif band == 'r':
            color = 'C1'
        elif band == 'i':
            color = 'C3'
        elif band == 'z':
            color = 'C5'
        else:
            color = 'k'  # should be y
        src = butler.get("src", dataId={"visit": visit, "detector": detectorToUse},
                         collections=collections, instrument=instrument)
        flag_src = [False]*len(src)
        f_nJy = src['base_PsfFlux_instFlux']
        for flag in badFlagList:
            flag_src |= src[flag]

        good_fluxes = np.array([s for s, f in zip(f_nJy, flag_src) if ~f and s > 0])
        mags = (good_fluxes*u.nJy).to_value(u.ABmag)

        if min_magnitude is None:
            min_magnitude = np.floor(np.min(mags)/binsize)*binsize
        if max_magnitude is None:
            max_magnitude = np.ceil(np.max(mags)/binsize)*binsize
        nbins = int((max_magnitude - min_magnitude)/binsize)
        hist, bin_edges = np.histogram(mags, bins=nbins, range=(min_magnitude, max_magnitude))
        bin_centers = [(bin_edges[i] + bin_edges[i + 1])/2 for i in range(len(bin_edges) - 1)]
        plt.plot(bin_centers, hist/area/binsize, label=visit, color=color)
    plt.title(f'{bandToPlot} Source Counts')
    plt.xlabel('Magnitude')
    plt.ylabel('Detected sources per mag per deg^2')
    ax.set_yscale('log')