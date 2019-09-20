import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import sqlite3
from astropy.visualization import (ZScaleInterval, SqrtStretch, ImageNormalize)

import lsst.daf.persistence as dafPersist
import lsst.afw.display as afwDisplay
import lsst.geom

"""Script to make light curve plots for DIAObjects using a Prompt Products
Database (PPDB) resulting from a run of ap_pipe.

Two plots are created: one for the light curve with a single set of
processed/template/difference cutouts, and one showing all of the
processed/difference cutouts.
"""


def main():
    afwDisplay.setDefaultBackend('matplotlib')
    script = sys.argv[0]
    try:
        repo = sys.argv[1]
        templateRepo = sys.argv[2]
    except:
        print('###')
        print('Run this script with 2-3 arguments: repo (required), \
               templateRepo (required), and dbName (optional)')
        print('Please note that dbName path must be relative to repo')
        print('For example, python {0} rerun/MyAmazingRerun rerun/TemplateLand association.db'.format(script))
        print('###')
    else:
        repo = os.path.normpath(repo)
        rerunName = os.path.basename(repo)
        # TODO: UN-HARDWIRE THIS
        miniRegion = 'and ra > 155.2 and ra < 155.3 and decl < -5.6 and decl > -5.8 \
                      and nDiaSources > 5'  # and flags == 0'
        try:
            dbName = sys.argv[3]
        except:
            print('Using default dbName, association.db')
            dbName = 'association.db'
        finally:
            print('Loading PPDB Objects...')
            objTable = loadPpdbObjects(repo, dbName, filter=miniRegion)
            objIdList = list(objTable['diaObjectId'])
            print('Loaded {0} DIA Objects'.format(len(objIdList)))
            # if len(objIdList) > 20:
            #     objIdList = objIdList[0:20]
            #     print('WARNING: objIdList has been truncated!')
            print('Plotting light curves for {0} objects...'.format(len(objIdList)))
            patchList = ['11,8', '11,9', '11,10']  # TODO: UN-HARDWIRE THIS
            pdfLightcurves = PdfPages(rerunName + '_lcs.pdf')
            pdfCutouts = PdfPages(rerunName + '_cos.pdf')
            for obj in objIdList:
                plotLightcurve(obj, objTable, repo, dbName, templateRepo,
                               patchList, plotAllCutouts=True,
                               pdfLc=pdfLightcurves, pdfCo=pdfCutouts)
            pdfLightcurves.close()
            pdfCutouts.close()


def in_ipynb():
    """Determine whether the code is being run in a notebook.
    """
    try:
        cfg = get_ipython().config
        if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
            return True
        elif 'LazyConfigValue' in str(cfg['IPKernelApp']['parent_appname']):
            return True
        else:
            return False
    except NameError:
        return False


def loadPpdbObjects(repo, dbName='association.db', filter=''):
    """Load select DIAObject columns from a PPDB into a pandas dataframe.

    Parameters
    ----------
    repo : `str`
        Path to an output repository from an ap_pipe run.
    dbName : `str`, optional
        Name of the PPDB, which must reside in (or relative to) repo.
    filter : `str`, optional
        Criteria defining which objects to load from the PPDB.
        The default is no additional filter (i.e., loading all objects);
        be careful - we will be plotting each object, so we want a number
        we can actually work with!

    Returns
    -------
    objTable : `pandas.DataFrame`
        DIA Object Table.
    """
    connection = sqlite3.connect(os.path.join(repo, dbName))

    # These are the tables available in the ppdb
    tables = {'obj': 'DiaObject', 'src': 'DiaSource', 'ccd': 'CcdVisit'}

    # Only get objects with validityEnd NULL because that means they are still valid
    objTable = pd.read_sql_query('select diaObjectId, ra, decl, nDiaSources, \
                                  gPSFluxMean, validityEnd, flags from {0} \
                                  where validityEnd is NULL \
                                  {1};'.format(tables['obj'], filter), connection)
    return objTable


def loadPpdbSources(dbPath, obj):
    """Load select DIAObject columns from a PPDB into a pandas dataframe.

    Parameters
    ----------
    repo : `str`
        Path to the PPDB.
    obj : `int`
        DIA Object for which we want to retrieve constituent DIA Sources.

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Object Table containing only objects with validityEnd NULL.
        Columns selected are presently hard-wired here.
    """
    connection = sqlite3.connect(dbPath)

    # These are the tables available in the ppdb
    tables = {'obj': 'DiaObject', 'src': 'DiaSource', 'ccd': 'CcdVisit'}

    # Load all information needed for light curves
    srcTable = pd.read_sql_query('select diaSourceId, diaObjectId, \
                                  ra, decl, ccdVisitId, \
                                  midPointTai, apFlux, psFlux, apFluxErr, \
                                  psFluxErr, totFlux, totFluxErr, flags from {0} \
                                  where diaObjectId = {1};'.format(tables['src'], obj), connection)
    return srcTable


def getTemplateCutout(scienceImage, templateRepo, centerSource, size=lsst.geom.Extent2I(30, 30),
                      templateDataType='deepCoadd', filter='g', templateVisit=None):
    """Retrieve cutout of difference imaging template.

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
    if 'Coadd' in templateDataType or 'coadd' in templateDataType:
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
                templateCutout = coaddTemplate.getCutout(centerSource, size)
            except:
                continue
            else:
                # Only loop through patches until you find one containing the source
                break
    else:  # it's a calexp or an instcal; assume DECam for now
        ccdnum = int(str(scienceImage.getInfo().getVisitInfo().getExposureId())[-2:])  # sorry
        templateDataId = dict(
            filter=filter,
            visit=templateVisit,
            ccdnum=ccdnum,
        )
        templateExposure = templateButler.get(templateDataType, dataId=templateDataId)
        templateCutout = templateExposure.getCutout(centerSource, size)

    print("Template dataId %s" % templateDataId)
    return templateCutout


def plotLightcurve(obj, objTable, repo, dbName, templateRepo,
                   useTotFlux=False, plotAllCutouts=False,
                   cutoutIdx=0, labelCutouts=False,
                   diffimType='deepDiff_differenceExp',
                   pdfLc=None, pdfCo=None, diffimRepo=None,
                   templateDataType='deepCoadd', templateVisitList=None,
                   orderVisits=True):
    """Plot lightcurve and processed, template, and difference image cutouts
    for one DIAObject. The lightcurve includes all associated DIASources.

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
        Location of PPDB relative to repo.
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
    print('Loading PPDB Sources...')
    dbPath = os.path.join(repo, dbName)
    srcTable = loadPpdbSources(dbPath, obj)
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
    except RuntimeError:
        calexpFirst = butler.get('instcal', dataIdDicts[cutoutIdx])
    calexpArray = calexpFirst.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
    calexpNorm = ImageNormalize(calexpArray, interval=ZScaleInterval(), stretch=SqrtStretch())
    plt.imshow(np.rot90(calexpArray), cmap='gray', norm=calexpNorm)

    # template image
    plt.subplot(232)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.title('Template', size=16)
    templateCutout = None
    if 'Coadd' in templateDataType or 'coadd' in templateDataType:
        templateCutout = getTemplateCutout(calexpFirst, templateRepo, centerSource)
        templateArray = templateCutout.getMaskedImage().getImage().getArray()
        templateNorm = ImageNormalize(templateArray, interval=ZScaleInterval(), stretch=SqrtStretch())
        plt.imshow(np.fliplr(templateArray), cmap='gray', norm=templateNorm)
    else:  # it's a calexp or an instcal, probably
        for visit in templateVisitList:
            try:
                templateCutout = getTemplateCutout(calexpFirst, templateRepo, centerSource,
                                                   templateDataType=templateDataType,
                                                   templateVisit=visit)
            except:
                continue  # loop through other possible visits
        templateArray = templateCutout.getMaskedImage().getImage().getArray()
        templateNorm = ImageNormalize(templateArray, interval=ZScaleInterval(), stretch=SqrtStretch())
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
    diffimNorm = ImageNormalize(diffimArray, interval=ZScaleInterval(), stretch=SqrtStretch())
    plt.imshow(np.rot90(diffimArray), cmap='gray', norm=diffimNorm)

    if plotAllCutouts:
        fig2 = plt.figure(figsize=(8, 8))  # optional figure with cutouts for all visits
        fig2.subplots_adjust(hspace=0, wspace=0)
        fig2.text(0.5, 0.9, str(obj), horizontalalignment='center', verticalalignment='top')
        for idx, dataId in enumerate(dataIdDicts):
            try:
                calexp = butler.get('calexp', dataId)
            except RuntimeError:
                calexp = butler.get('instcal', dataId)
            calexpArray = calexp.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
            calexpNorm = ImageNormalize(calexpArray, interval=ZScaleInterval(), stretch=SqrtStretch())
            if diffimRepo is not None:
                butlerDiffim = dafPersist.Butler(diffimRepo)
                diffim = butlerDiffim.get(diffimType, dataId)
            else:
                diffim = butler.get(diffimType, dataId)
            diffimArray = diffim.getCutout(centerSource, size).getMaskedImage().getImage().getArray()
            diffimNorm = ImageNormalize(diffimArray, interval=ZScaleInterval(), stretch=SqrtStretch())
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
        if pdfCo and not in_ipynb():
            fig2.savefig(pdfCo, format='pdf')
        elif not in_ipynb():
            fig2.savefig(str(obj) + '_co.png')
        else:
            plt.show()
    if pdfLc and not in_ipynb():
        fig1.savefig(pdfLc, format='pdf')
    elif not in_ipynb():
        fig1.savefig(str(obj) + '_lc.png')
    else:
        plt.show()


if __name__ == '__main__':
    main()
