import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import astropy.visualization as aviz

import lsst.daf.persistence as dafPersist
import lsst.geom
import lsst.pex.exceptions

import diaObjectAnalysis as doa

"""Utilities to make cutouts and light curves for DIA Sources
and DIA Objects using an Alert Production Database (APDB).
"""


def loadSelectApdbSources(dbName, diaObjectId, dbType='sqlite', schema=None):
    """Load select columns from DIASources for a single DIAObject
    from an APDB into a pandas dataframe.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, full filepath to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    diaObjectId : `int`
        DIA Object for which we want to retrieve constituent DIA Sources.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Source Table including the columns hard-wired below.
    """
    connection, tables = doa.connectToApdb(dbName, dbType, schema)

    # Load data from the source table
    srcTable = pd.read_sql_query('select "diaSourceId", "diaObjectId", \
                                  "ra", "decl", "ccdVisitId", "filterName", \
                                  "midPointTai", "apFlux", "psFlux", "apFluxErr", \
                                  "psFluxErr", "totFlux", "totFluxErr", "x", "y", \
                                  "ixxPSF", "iyyPSF", "ixyPSF", "flags" from {0}; \
                                  where "diaObjectId" = {1};'.format(tables['src'], diaObjectId), connection)
    return srcTable


def loadExposures(butler, dataId, collections, diffName='deep'):
    """Load a science exposure, difference image, and warped template.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.
    dataId : `dict`-like
        Gen3 data ID specifying at least instrument, visit, and detector.
    collections : `str` or `list`
        Gen3 collection or collections from which to load the exposures.
    diffName : `str`, optional
        Default is 'deep', but 'goodSeeing' may be needed instead.

    Returns
    -------
    science : `lsst.afw.Exposure`
        calexp corresponding to dataId and collections.
    difference : `lsst.afw.Exposure`
        differenceExp corresponding to dataId and collections.
    template : `lsst.afw.Exposure`
        warpedExp corresponding to dataId and collections.
    """
    science = butler.get('calexp', dataId=dataId, collections=collections)
    difference = butler.get(diffName + 'Diff_differenceExp', dataId=dataId, collections=collections)
    template = butler.get(diffName + 'Diff_warpedExp', dataId=dataId, collections=collections)
    return science, difference, template


def retrieveCutouts(butler, dataId, collections, center, size=lsst.geom.Extent2I(30, 30), diffName='deep'):
    """Return small cutout exposures for a science exposure, difference image,
    and warped template.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler in the repository corresponding to the output of an ap_pipe run.
    dataId : `dict`-like
        Gen3 data ID specifying at least instrument, visit, and detector.
    collections : `str` or `list`
        Gen3 collection or collections from which to load the exposures.
    center : `lsst.geom.SpherePoint`
        Desired center coordinate of cutout.
    size : `lsst.geom.Extent`, optional
        Desired size of cutout, default is 30x30 pixels
    diffName : `str`, optional
        Default is 'deep', but 'goodSeeing' may be needed instead.

    Returns
    -------
    scienceCutout : `lsst.afw.Exposure`
        Cutout of calexp at location 'center' of size 'size'.
    differenceCutout : `lsst.afw.Exposure`
        Cutout of deepDiff_differenceExp at location 'center' of size 'size'.
    templateCutout : `lsst.afw.Exposure`
        Cutout of deepDiff_warpedExp at location 'center' of size 'size'.
    """
    science, difference, template = loadExposures(butler, dataId, collections, diffName)
    scienceCutout = science.getCutout(center, size)
    differenceCutout = difference.getCutout(center, size)
    templateCutout = template.getCutout(center, size)
    return scienceCutout, differenceCutout, templateCutout


def plotCutout(scienceCutout, differenceCutout, templateCutout, output=None):
    """Plot the cutouts for one DIASource in one image.

    Parameters
    ----------
    scienceCutout : `lsst.afw.Exposure`
        Cutout of calexp returned by retrieveCutouts.
    differenceCutout : `lsst.afw.Exposure`
        Cutout of deepDiff_differenceExp returned by retrieveCutouts.
    templateCutout : `lsst.afw.Exposure`
        Cutout of deepDiff_warpedExp returned by retrieveCutouts.
    output : `str`, optional
        If provided, save png to disk at output filepath.
    """
    def do_one(ax, data, name):
        interval = aviz.ZScaleInterval()
        if name == 'Difference':
            norm = aviz.ImageNormalize(data, stretch=aviz.LinearStretch())
        else:
            norm = aviz.ImageNormalize(data, interval=interval, stretch=aviz.AsinhStretch(a=0.01))
        ax.imshow(data, cmap=cm.bone, interpolation="none", norm=norm)
        ax.axis('off')
        ax.set_title(name)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    do_one(ax1, templateCutout.image.array, "Template")
    do_one(ax2, scienceCutout.image.array, "Science")
    do_one(ax3, differenceCutout.image.array, "Difference")
    plt.tight_layout()

    if output is not None:
        plt.savefig(output, bbox_inches="tight")
        plt.close()


def getTemplateCutoutGen2(scienceImage, templateRepo, centerSource, size=lsst.geom.Extent2I(30, 30),
                          templateDataType='deepCoadd', filter='g', templateVisit=None):
    """Retrieve cutout of difference imaging template.

    This function is not maintained and is for use with gen2 middleware only.

    Parameters
    ----------
    scienceImage : `lsst.afw.image`
        Image which has been differenced that you want to get the template for.
    templateRepo : `str`
        Path to repository containing difference imaging templates.
    centerSource : `lsst.geom.SpherePoint`
        Desired center coordinate of cutout.
    size : `lsst.geom.Extent`, optional
        Desired size of cutout, default is 30x30 pixels
    templateDataType : `str`, optional
        Default is 'deepCoadd'. Other possibilities could be 'calexp', 'instcal', etc.
    filter : `str`, optional
        Filter to use in constructing template dataId. Default is 'g'.
    templateVisit : `int` or None, optional
        Must specify templateVisit if templateDataType is not some kind of coadd.

    Returns
    -------
    coaddCutout : `lsst.afw.image`
        A small image postage stamp cutout around the template source.
    """
    templateButler = dafPersist.Butler(templateRepo)
    skyMap = templateButler.get(datasetType=templateDataType + '_skyMap')
    expWcs = scienceImage.getWcs()
    expBoxD = lsst.geom.Box2D(scienceImage.getBBox())
    expBoxD.grow(10)
    ctrSkyPos = expWcs.pixelToSky(expBoxD.getCenter())
    tractInfo = skyMap.findTract(ctrSkyPos)
    skyCorners = [expWcs.pixelToSky(pixPos) for pixPos in expBoxD.getCorners()]
    patchList = tractInfo.findPatchList(skyCorners)
    for patchInfo in patchList:
        templateDataId = dict(
            tract=tractInfo.getId(),
            patch="%s,%s" % (patchInfo.getIndex()[0], patchInfo.getIndex()[1]),
            filter=filter,
        )
        try:
            coaddTemplate = templateButler.get('deepCoadd', dataId=templateDataId)
        except dafPersist.butlerExceptions.NoResults:
            template = None
            continue
        else:
            try:
                template = coaddTemplate.getCutout(centerSource, size)
            except lsst.pex.exceptions.InvalidParameterError:
                template = None
                continue
            # Only loop through patches until you find one containing the source
            break

    print("Template dataId %s" % templateDataId)
    return template


def plotLightcurveGen2(obj, objTable, repo, dbName, templateRepo,
                       useTotFlux=False, plotAllCutouts=False,
                       cutoutIdx=0, labelCutouts=False,
                       diffimType='deepDiff_differenceExp',
                       pdfLc=None, pdfCo=None, diffimRepo=None,
                       templateDataType='deepCoadd', templateVisitList=None,
                       orderVisits=True):
    """Plot lightcurve and processed, template, and difference image cutouts
    for one DIAObject. The lightcurve includes all associated DIASources.

    This function is not maintained and is for use with gen2 middleware only.

    Nothing is returned, but a plot is saved to disk.

    Parameters
    ----------
    obj : `int`
        DIA Object ID.
    objTable : `pandas.DataFrame`
        DIA Object Table which contains obj.
    repo : `str`
        Path to repository containing processed visit images.
    dbName : `str`
        Location of APDB relative to repo.
    templateRepo : `str`
        Path to repository containing templates for image differencing.
    useTotFlux : `bool`, optional, default False
    plotAllCutouts : `bool`, optional, default False
    cutoutIdx : `int`, optional, default 0
    labelCutouts : `bool`, optional, default False
    diffimType : `str`, optional, default deep (not dcr)
        Typically either 'deepDiff_differenceExp' or 'dcrDiff_differenceExp'
    pdfLc : `str`
        Filename to save PDF light curve plots to
    pdfCo : `str`
        Filename to save PDF cutout plots to
    diffimRepo : `str`, optional
        Path to repository containing diffims. Only use if this differs from repo.
    templateDataType : `str`, optional
        Specify, e.g., 'calexp' or 'instcal' if you don't have a 'deepCoadd'.
    templateVisitList : `list`, optional
        If you're using non-coadd templates and want to retrieve a
        cutout, provide a list of possible template visits to try here.
    """
    print('Loading APDB Sources...')
    dbPath = os.path.join(repo, dbName)
    srcTable = loadSelectApdbSources(dbPath, obj)
    if orderVisits:
        srcTable = srcTable.sort_values("ccdVisitId")
    ra = objTable.loc[objTable['diaObjectId'] == obj, 'ra']
    dec = objTable.loc[objTable['diaObjectId'] == obj, 'decl']
    #  flags = srcTable['flags']
    dataIds = srcTable['ccdVisitId'].values  # these are ints
    dataIdDicts = []
    for dataId in dataIds:
        visit = int(str(dataId)[0:6])
        ccdnum = int(str(dataId)[6:])
        dataIdDict = {'visit': visit, 'ccdnum': ccdnum}
        dataIdDicts.append(dataIdDict)
    centerSource = lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees)
    size = lsst.geom.Extent2I(30, 30)

    print('DIAObject ID:', obj)
    # print('Flags:', flags)
    print('RA (deg):', ra.values)
    print('Dec (deg):', dec.values)
    print('Number of DIASources:', len(srcTable['diaSourceId'].values))
    # print('DIASource IDs:', srcTable['diaSourceId'].values)
    # print('Data IDs:', dataIdDicts)

    fig1 = plt.figure()
    fig1.text(0.5, 0.5, str(obj), horizontalalignment='center', verticalalignment='center')

    # light curve with psFlux by default (uses totFlux if useTotFlux=True)
    plt.subplot(212)
    plt.xlabel('Time (MJD)', size=16)
    if not useTotFlux:
        plt.errorbar(srcTable['midPointTai'], srcTable['psFlux']*1e9, yerr=srcTable['psFluxErr']*1e9,
                     ls=':', marker='o', color='#2979C1')
        plt.ylabel('Difference Flux (nJy)', size=16)
    else:
        plt.errorbar(srcTable['midPointTai'], srcTable['totFlux']*1e9, yerr=srcTable['totFluxErr']*1e9,
                     ls=':', marker='o', color='#2979C1')
        plt.ylabel('Flux (nJy)', size=16)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # processed image
    plt.subplot(231)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.title('Processed', size=16)
    butler = dafPersist.Butler(repo)
    try:
        calexpFirst = butler.get('calexp', dataIdDicts[cutoutIdx])
    except dafPersist.butlerExceptions.NoResults:
        calexpFirst = butler.get('instcal', dataIdDicts[cutoutIdx])
    calexpArray = calexpFirst.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
    calexpNorm = aviz.ImageNormalize(calexpArray,
                                     interval=aviz.ZScaleInterval(),
                                     stretch=aviz.SqrtStretch())
    plt.imshow(np.rot90(calexpArray), cmap='gray', norm=calexpNorm)

    # template image
    plt.subplot(232)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.title('Template', size=16)
    templateCutout = None
    if 'Coadd' in templateDataType or 'coadd' in templateDataType:
        templateCutout = getTemplateCutoutGen2(calexpFirst, templateRepo, centerSource)
        templateArray = templateCutout.getMaskedImage().getImage().getArray()
        templateNorm = aviz.ImageNormalize(templateArray,
                                           interval=aviz.ZScaleInterval(),
                                           stretch=aviz.SqrtStretch())
        plt.imshow(np.fliplr(templateArray), cmap='gray', norm=templateNorm)
    else:  # it's a calexp or an instcal, probably
        for visit in templateVisitList:
            try:
                templateCutout = getTemplateCutoutGen2(calexpFirst, templateRepo, centerSource,
                                                       templateDataType=templateDataType,
                                                       templateVisit=visit)
            except dafPersist.butlerExceptions.NoResults:
                continue  # loop through other possible visits
        templateArray = templateCutout.getMaskedImage().getImage().getArray()
        templateNorm = aviz.ImageNormalize(templateArray,
                                           interval=aviz.ZScaleInterval(),
                                           stretch=aviz.SqrtStretch())
        plt.imshow(np.rot90(templateArray), cmap='gray', norm=templateNorm)

    # difference image
    plt.subplot(233)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.title('Difference', size=16)
    if diffimRepo is not None:
        butlerDiffim = dafPersist.Butler(diffimRepo)
        diffimFirst = butlerDiffim.get(diffimType, dataIdDicts[cutoutIdx])
    else:
        diffimFirst = butler.get(diffimType, dataIdDicts[cutoutIdx])
    diffimArray = diffimFirst.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
    diffimNorm = aviz.ImageNormalize(diffimArray,
                                     interval=aviz.ZScaleInterval(),
                                     stretch=aviz.SqrtStretch())
    plt.imshow(np.rot90(diffimArray), cmap='gray', norm=diffimNorm)

    if plotAllCutouts:
        fig2 = plt.figure(figsize=(8, 8))  # optional figure with cutouts for all visits
        fig2.subplots_adjust(hspace=0, wspace=0)
        fig2.text(0.5, 0.9, str(obj), horizontalalignment='center', verticalalignment='top')
        for idx, dataId in enumerate(dataIdDicts):
            try:
                calexp = butler.get('calexp', dataId)
            except dafPersist.butlerExceptions.NoResults:
                calexp = butler.get('instcal', dataId)
            calexpArray = calexp.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
            calexpNorm = aviz.ImageNormalize(calexpArray,
                                             interval=aviz.ZScaleInterval(),
                                             stretch=aviz.SqrtStretch())
            if diffimRepo is not None:
                butlerDiffim = dafPersist.Butler(diffimRepo)
                diffim = butlerDiffim.get(diffimType, dataId)
            else:
                diffim = butler.get(diffimType, dataId)
            diffimArray = diffim.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
            diffimNorm = aviz.ImageNormalize(diffimArray,
                                             interval=aviz.ZScaleInterval(),
                                             stretch=aviz.SqrtStretch())
            plt.subplot(12, 12, idx+1)
            plt.gca().get_xaxis().set_ticks([])
            plt.gca().get_yaxis().set_ticks([])
            plt.imshow(np.rot90(calexpArray), cmap='gray', norm=calexpNorm)
            if labelCutouts:
                if idx == 0:
                    plt.text(1, 26, 'Proc', color='lime', size=8)
                plt.text(2, 5, str(srcTable['midPointTai'][idx])[1:8], color='lime', size=8)
                # TODO: PRINT FLAG INFORMATION ON/NEAR CUTOUTS...?
            plt.subplot(12, 12, idx+72+1)
            plt.gca().get_xaxis().set_ticks([])
            plt.gca().get_yaxis().set_ticks([])
            plt.imshow(np.rot90(diffimArray), cmap='gray', norm=diffimNorm)
            if labelCutouts:
                if idx == 0:
                    plt.text(1, 26, 'Diff', color='lime', size=8)
                plt.text(2, 5, str(srcTable['midPointTai'][idx])[1:8], color='lime', size=8)
        if pdfCo:  # may not work in notebook environment
            fig2.savefig(pdfCo, format='pdf')
        else:
            plt.show()
    if pdfLc:  # may not work in notebook environment
        fig1.savefig(pdfLc, format='pdf')
    else:
        plt.show()
