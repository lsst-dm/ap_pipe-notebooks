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
import lsst.afw.display as afwDisplay
import lsst.geom
from lsst.ap.association import MapDiaSourceConfig, UnpackApdbFlags
import lsst.afw.cameraGeom as cameraGeom

from diaObjectAnalysis import loadAllApdbObjects, loadAllApdbSources
from plotLightcurve import getTemplateCutout
"""
Collection of plots that can be made using info in the APDB.

In some cases, additional inputs may be needed, and that will be
clearly indicated in the docstrings.

Future plots could include:
- plot DIA Sources on focal plane on one visit or ccd at a time
- plot various totFlux vs psFlux and psFlux vs apFlux (or magnitudes)
- plot stats about the images comprising the template
"""


def addVisitCcdToSrcTable(sourceTable, cam='decam'):
    """Add visit and ccd columns to sourceTable dataframe.
    This currently works for DECam and HSC only.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    cam : `str`, one of either 'decam' or 'hsc', default is 'decam'
        Needed to properly add the "ccd" and "visit" columns to the sourceTable.

    Returns
    -------
    sourceTable : `pandas.core.frame.DataFrame`
        The same as the input sourceTable, with new visit and ccd columns.
    """
    if cam == 'decam':
        sourceTable['ccd'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[-2:])
        sourceTable['visit'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[:-2])
    elif cam == 'hsc':
        sourceTable['visit'] = sourceTable.ccdVisitId.apply(lambda x: int(x/200.))
        sourceTable['ccd'] = sourceTable.ccdVisitId.apply(lambda x: int(np.round((x/200. - int(x/200.))*200)))
    else:
        raise ValueError('This plotting utility only works for DECam and HSC.')
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


def plotDiaSourcesPerVisit(repo, sourceTable, title=''):
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
    """
    ccdArea, visitArea = getCcdAndVisitSizeOnSky(repo, sourceTable)
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
    title : `str`
        Title for the plot, optional.
    """
    sourceTable['date_time'] = pd.to_datetime(sourceTable.midPointTai,
                                              unit='D',
                                              origin=pd.Timestamp('1858-11-17'))
    sourceTable['date'] = sourceTable['date_time'].dt.date
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


def getCcdCorners(butler, sourceTable):
    """Get corner coordinates for a range of ccds.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.

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
            bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
        except dafPersist.NoResults:  # silently skip any ccds that don't exist
            for visit in visits[1:]:  # try the other visits just in case
                visit = int(visit)
                try:
                    bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
                except dafPersist.NoResults:
                    continue
                break
        else:
            camera = butler.get('camera')
            cornerList.append([visit, ccd] + [val for pair in bbox.getCorners()
                              for val in ccd2focalPlane(pair[0], pair[1], ccd, camera)])
    corners = pd.DataFrame(cornerList)
    corners['width'] = corners[6] - corners[8]
    corners['height'] = corners[7] - corners[5]
    return corners


def getCcdAndVisitSizeOnSky(repo, sourceTable):
    butler = dafPersist.Butler(repo)
    visits = np.unique(sourceTable.visit)
    ccds = np.unique(sourceTable.ccd)
    nGoodCcds = len(ccds)
    calexp = butler.get('calexp', visit=int(visits[0]), ccd=int(ccds[0]))
    bbox = butler.get('calexp_bbox', visit=int(visits[0]), ccd=int(ccds[0]))
    pixelScale = calexp.getWcs().getPixelScale().asArcseconds()
    ccdArea = (pixelScale*pixelScale*bbox.getArea()*u.arcsec**2).to(u.deg**2).value
    visitArea = ccdArea * nGoodCcds
    return ccdArea, visitArea


def plotDiaSourceDensityInFocalPlane(repo, sourceTable, cmap=mpl.cm.Blues, title=''):
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
    """
    ccdArea, visitArea = getCcdAndVisitSizeOnSky(repo, sourceTable)
    ccds = np.unique(sourceTable['ccd'])
    ccdMax = int(np.max(ccds))
    nVisits = len(np.unique(sourceTable['visit'].values))
    ccdGroup = sourceTable.groupby('ccd')
    ccdSourceCount = ccdGroup.visit.count().values/nVisits/ccdArea
    # DIA Source count per visit per square degree, for each CCD
    butler = dafPersist.Butler(repo)
    corners = getCcdCorners(butler, sourceTable)
    norm = mpl.colors.Normalize(vmin=np.min(ccdSourceCount), vmax=np.max(ccdSourceCount))

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111, aspect='equal')

    for index, row in corners.iterrows():
        try:
            if ccdMax < 100:  # it's DECam
                averageFocalPlane = ccdGroup.get_group('%02d' % row[1]).x.count()/nVisits/ccdArea
            else:  # it's HSC
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
               cam='decam'):
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
    badFlagList :  `list`
        Names of flags presumed to each indicate a DIA Source is garbage.
    cam : `str`, one of either 'decam' or 'hsc', default is 'decam'
        Needed to properly add the "ccd" and "visit" columns to the sourceTable.

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

    objTable = loadAllApdbObjects(dbPath, dbType=dbType)
    srcTable = loadAllApdbSources(dbPath, dbType=dbType)
    srcTable = addVisitCcdToSrcTable(srcTable, cam=cam)
    flagTable, srcTableFlags, flagFilter, \
        goodSrc, goodObj = makeSrcTableFlags(srcTable, objTable, badFlagList=badFlagList)
    return objTable, srcTable, goodObj, goodSrc


def makeSrcTableFlags(sourceTable, objectTable,
                      badFlagList=['base_PixelFlags_flag_bad',
                                   'base_PixelFlags_flag_suspect',
                                   'base_PixelFlags_flag_saturatedCenter']):
    """Apply flag filters to a DIA Source and a DIA Object table.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    objectTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Objects from the same APDB.
    badFlagList :  `list`
        Names of flags presumed to each indicate a DIA Source is garbage.

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
    config = MapDiaSourceConfig()
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


def plotSeeingHistogram(repo, sourceTable, ccd=35):
    """Plot distribution of visit seeing.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    ccd : `int`
        The ccd being considered, default 35.
    """
    fwhm = pd.DataFrame()
    butler = dafPersist.Butler(repo)
    visits = np.unique(sourceTable['visit'])
    radii = []
    for visit in visits:
        miniCalexp = butler.get('calexp_sub', visit=int(visit),
                                ccd=ccd, bbox=lsst.geom.Box2I())
        psf = miniCalexp.getPsf()
        psfSize = psf.computeShape().getDeterminantRadius()
        radii.append(psfSize*2.355)  # convert sigma to FWHM
    fwhm['visit'] = pd.Series(visits)
    fwhm['radius'] = pd.Series(radii, index=fwhm.index)
    pixelScale = miniCalexp.getWcs().getPixelScale().asArcseconds()  # same for all visits
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.hist(fwhm['radius'].values, alpha=0.5)
    plt.xlabel('Seeing FWHM (pixels)')
    plt.ylabel('Visit count')
    secax = ax.secondary_xaxis('top', functions=(lambda x: x*pixelScale, lambda x: x/pixelScale))
    secax.set_xlabel('Seeing FWHM (arcseconds)')


def plotDiaSourcesInFocalPlane(repo, sourceTable, gridsize=(400, 400), title=''):
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
    """
    butler = dafPersist.Butler(repo)
    camera = butler.get('camera')
    corners = getCcdCorners(butler, sourceTable)
    xFP_list = []
    yFP_list = []
    for index, row in sourceTable.iterrows():
        xFP, yFP = ccd2focalPlane(row['x'], row['y'], row['ccd'], camera=camera)
        xFP_list.append(xFP)
        yFP_list.append(yFP)
    sourceTable['xFP'] = pd.Series(xFP_list, index=sourceTable.index)
    sourceTable['yFP'] = pd.Series(yFP_list, index=sourceTable.index)

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111, aspect='equal')

    for index, row in corners.iterrows():
        ax1.add_patch(patches.Rectangle((row[7], row[6]),
                                        -row.height,
                                        -row.width,
                                        fill=False))
        ax1.text(row[7] - row.height/2, row[6] - row.width/2, '%d' % (row[1]))
        plt.plot(row[7] - row.height/2, row[6] - row.width/2, ',')

    ax1.hexbin(sourceTable['yFP'], sourceTable['xFP'], gridsize=gridsize, bins='log', cmap='Blues')
    ax1.set_title('DIA Sources in focal plane coordinates %s' % (title))

    ax1.set_xlabel('Focal Plane X', size=16)
    ax1.set_ylabel('Focal Plane Y', size=16)
    ax1.invert_yaxis()
    ax1.invert_xaxis()


def plotDiaSourcesOnSkyGrid(repo, sourceTable, title=None, color='C0'):
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
                   marker='.', s=0.1, alpha=0.2)
        ax.set_title(visit, size=8)
        ax.invert_xaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0)
    if title:
        fig.suptitle(title)
