import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import functools
import operator

import lsst.daf.persistence as dafPersist
import lsst.afw.display as afwDisplay
import lsst.geom
from lsst.ap.association import MapDiaSourceConfig, UnpackApdbFlags
import lsst.afw.cameraGeom as cameraGeom

from diaObjectAnalysis import loadAllApdbObjects, loadAllApdbSources
import plotLightcurve as plc
"""
Collection of plots that can be made using info in the APDB.

In some cases, additional inputs may be needed, and that will be
clearly indicated in the docstrings.

Future plots could include:
- plot DIA Sources on focal plane on one visit or ccd at a time
- plot various totFlux vs psFlux and psFlux vs apFlux (or magnitudes)
- plot stats about the images comprising the template
"""


def addVisitCcdToSrcTable(sourceTable):
    """Add visit and ccd columns to sourceTable dataframe.

    This assumes the `ccdVisitId` field may be cast to type `str` and
    is arranged such that the last 2 characters represent a ccd and all
    preceding characters represent a visit.

    TODO: generalize this for other exposure ID formats.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.

    Returns
    -------
    sourceTable : `pandas.core.frame.DataFrame`
        The same as the input sourceTable, with new visit and ccd columns.
    """
    sourceTable['ccd'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[-2:])
    sourceTable['visit'] = sourceTable.ccdVisitId.apply(lambda x: str(x)[:-2])
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


def plotDiaSourcesPerVisit(sourceTable, approxVisitArea=0.045, title=''):
    """Plot DIA Sources per visit.

    The plot will have two y-axes: number of DIA Sources per square degree and
    median FWHM per ixx or iyy in pixels.

    Parameters
    ----------
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    approxVisitArea : `float`
        Approximate area of one visit (all ccds), in square degrees.
        Default set to an appropriate DECam value.
    title : `str`
        Title for the plot, optional.
    """
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    traceRadius = np.sqrt(0.5) * np.sqrt(sourceTable.ixxPSF + sourceTable.iyyPSF)
    sourceTable['seeing'] = 2*np.sqrt(2*np.log(2)) * traceRadius
    visitGroup = sourceTable.groupby('visit')
    plot2axes(visitGroup.visit.first().values,
              visitGroup.ccd.count().values/approxVisitArea,
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


def ccd2focalPlane(x, y, ccd, butler):
    """Retrieve focal plane coordinates.

    Parameters
    ----------
    x : `int` or `float`
        X-coordinate from ccd bbox.
    y : `int` or `float`
        Y-coordinate from ccd bbox.
    ccd : `int`, or can be cast as int
        The ccd being considered.
    butler : `lsst.daf.persistence.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.

    Returns
    -------
    point[0] : `int` or `float`
        X-coordinate in focal plane.
    point[1] : `int` or `float`
        Y-coordinate in focal plane.
    """
    camera = butler.get('camera')
    detector = camera[int(ccd)]
    point = detector.transform(lsst.geom.Point2D(x, y),
                               cameraGeom.PIXELS, cameraGeom.FOCAL_PLANE)
    return point[0], point[1]


def getCcdCorners(butler, sourceTable, ccdMin=1, ccdMax=63):
    """Get corner coordinates for a range of ccds.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    ccdMin : `int`
        Starting value for ccd range, default 1 for DECam.
    ccdMax : `int`
        Ending value for ccd range (final ccd + 1), default 63 for DECam.

    Returns
    -------
    corners : `pandas.core.frame.DataFrame`
        Dataframe containing focal plane coordinates for all the ccd corners.
    """
    cornerList = []
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    visits = np.unique(sourceTable['visit'])
    visit = int(visits[0])  # shouldn't matter what visit we use
    for ccd in range(ccdMin, ccdMax):
        try:
            bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
        except dafPersist.NoResults:  # silently skip any ccds that don't exist
            # print('No bbox for visit, ccd ', visit, ccd)
            visit = int(visits[1])  # try a different visit just in case
            try:
                bbox = butler.get('calexp_bbox', visit=visit, ccd=ccd)
            except dafPersist.NoResults:
                # print('Also no bbox for visit, ccd ', visit, ccd)
                continue  # give up
            continue
        else:
            cornerList.append([visit, ccd] + [val for pair in bbox.getCorners()
                              for val in ccd2focalPlane(pair[0], pair[1], ccd, butler)])
    corners = pd.DataFrame(cornerList)
    corners['width'] = corners[6] - corners[8]
    corners['height'] = corners[7] - corners[5]
    return corners


def plotDiaSourceDensityInFocalPlane(repo, sourceTable, approxCcdArea=2.678,
                                     cmap=mpl.cm.Blues, title=''):
    """Plot average density of DIA Sources in the focal plane (per CCD).

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    approxCcdArea : `float`
        Approximate area per ccd, in square degrees.
        Default set to an appropriate DECam value.
    cmap : `matplotlib.colors.ListedColormap`
        Matplotlib colormap.
    title : `str`
        String to append to the plot title, optional.
    """
    butler = dafPersist.Butler(repo)
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    nVisits = len(np.unique(sourceTable['visit'].values))
    grp = sourceTable.groupby('ccd')
    yy = grp.visit.count().values/nVisits/approxCcdArea
    corners = getCcdCorners(butler, sourceTable)
    norm = mpl.colors.Normalize(vmin=np.min(yy), vmax=np.max(yy))

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111, aspect='equal')

    for index, row in corners.iterrows():
        color = cmap(norm(grp.get_group('%02d' % row[1]).x.count()/nVisits/approxCcdArea))
        ax1.add_patch(patches.Rectangle((row[7], row[6]),
                      -row.height,
                      -row.width,
                      fill=True,
                      color=color))
        ax1.text(row[7]-row.height/2, row[6]-row.width/2, '%d' % (row[1]), fontsize=14)
        plt.plot(row[7]-row.height/2, row[6]-row.width/2, ',')
    ax1.set_title('Mean DIA Source density in focal plane coordinates %s' % (title))
    ax1.set_xlabel('Focal Plane X')
    ax1.set_ylabel('Focal Plane Y')
    ax1 = plt.gca()
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, fraction=0.04, pad=0.04)
    cb.set_label('DIA Sources per sq. deg.', rotation=90)


def loadTables(repo, dbName='association.db', isVerify=False,
               badFlagList=['base_PixelFlags_flag_bad',
                            'base_PixelFlags_flag_suspect',
                            'base_PixelFlags_flag_saturatedCenter']):
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
    if not isVerify:  # APDB is in repo (default)
        dbPath = os.path.abspath(os.path.join(repo, dbName))
    else:  # APDB is one level above repo
        repoUpOne = os.path.dirname(repo)
        dbPath = os.path.abspath(os.path.join(repoUpOne, dbName))
    objTable = loadAllApdbObjects(dbPath)
    srcTable = loadAllApdbSources(dbPath)
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
    srcCcdVisitIds = sourceTable.loc[sourceTable['diaObjectId'] == objId, 'ccdVisitId'].values
    dataIds = [{'visit': int(str(dataId)[0:6]), 'ccdnum': int(str(dataId)[6:])} for dataId in srcCcdVisitIds]

    objLoc = lsst.geom.SpherePoint(objRa, objDec, lsst.geom.degrees)
    calexpFirst = butler.get('calexp', dataId=dataIds[0])
    templateCutout = plc.getTemplateCutout(calexpFirst, repo, objLoc)
    wcs = templateCutout.getWcs()

    disp = afwDisplay.Display(plotIdx, reopenPlot=True)
    disp.setMaskTransparency(100)
    disp.scale('asinh', 'zscale', Q=4)
    disp.mtv(templateCutout, title=objId)

    with disp.Buffering():  # obviously
        coordObj = wcs.skyToPixel(objLoc)
        disp.dot('+', *coordObj, ctype='C0', size=5)
        for dataId, srcId, ra, dec in zip(dataIds, srcIds, srcRas, srcDecs):
            calexp = butler.get('calexp', dataId=dataId)
            psf = calexp.getPsf()
            psfSize = psf.computeShape().getDeterminantRadius()
            coordSrc = wcs.skyToPixel(lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees))
            if goodSrc is not None:
                if srcId in goodSrc['diaSourceId'].values:
                    disp.dot('o', *coordSrc, ctype='C3', size=psfSize)  # green
                else:
                    disp.dot('o', *coordSrc, ctype='C2', size=psfSize)  # red
            else:
                disp.dot('o', *coordSrc, ctype='C0', size=psfSize)  # blue


def plotSeeingHistogram(repo, sourceTable, ccd=35, ccdXcenter=1000,
                        ccdYcenter=2000, pixelScale=0.26, title=''):
    """Plot distribution of visit seeing.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    ccd : `int`
        The ccd being considered.
    ccdXcenter : `int`
        Pixel center of ccd in x-direction (e.g., 1000 for 2000x4000 ccd).
    ccdYcenter : `int`
        Pixel center of ccd in y-direction (e.g., 2000 for 2000x4000 ccd).
    pixelScale : `float`
        Pixel scale of the instrument in arcseconds per pixel, default 0.26.
    title : `str`
        Title for the plot, optional.
    """
    fwhm = pd.DataFrame()
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    butler = dafPersist.Butler(repo)
    visits = np.unique(sourceTable['visit'])
    radii = []
    for visit in visits:
        miniCalexp = butler.get('calexp_sub', visit=int(visit),
                                ccd=ccd, bbox=lsst.geom.Box2I())
        psf = miniCalexp.getPsf()
        shape = psf.computeShape(lsst.geom.Point2D(ccdXcenter, ccdYcenter))
        radii.append(shape.getDeterminantRadius()*2.355)  # convert sigma to FWHM
    fwhm['visit'] = pd.Series(visits)
    fwhm['radius'] = pd.Series(radii, index=fwhm.index)

    plt.figure(figsize=(6, 8))
    plt.subplot(211)
    plt.hist(fwhm['radius'].values*pixelScale, alpha=0.5)
    plt.xlabel('Seeing FWHM (arcseconds)')
    plt.ylabel('Visit count')
    plt.subplot(212)
    plt.hist(fwhm['radius'].values, alpha=0.5)
    plt.xlabel('Seeing FWHM (pixels)')
    plt.ylabel('Visit count')


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
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    butler = dafPersist.Butler(repo)
    corners = getCcdCorners(butler, sourceTable, ccdMin=1, ccdMax=63)
    xFP_list = []
    yFP_list = []
    for index, row in sourceTable.iterrows():
        xFP, yFP = ccd2focalPlane(row['x'], row['y'], row['ccd'], butler=butler)
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

    ax1.set_xlabel('Focal Plane X')
    ax1.set_ylabel('Focal Plane Y')
    ax1.invert_yaxis()
    ax1.invert_xaxis()


def plotDiaSourcesOnSkyGrid(repo, sourceTable):
    """Make a multi-panel plot of DIA Sources for each visit on the sky.

    Parameters
    ----------
    repo : `str`
        Repository corresponding to the output of an ap_pipe run.
    sourceTable : `pandas.core.frame.DataFrame`
        Pandas dataframe with DIA Sources from an APDB.
    title : `str`
        String to append to the plot title, optional.
    """
    sourceTable = addVisitCcdToSrcTable(sourceTable)
    fig = plt.figure(figsize=(9, 10))
    for count, visit in enumerate(np.unique(sourceTable['visit'].values)):
        idx = sourceTable.visit == visit
        ax = fig.add_subplot(10, 10, count + 1, aspect='equal')
        ax.scatter(sourceTable.ra[idx], sourceTable.decl[idx], c='C0',
                   marker='.', s=0.1, alpha=0.2)
        ax.set_title(visit, size=8)
        ax.invert_xaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0)
