import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
import sqlite3
import psycopg2

"""Script to make plots of DIAObjects using a Alert Production Database
(APDB) resulting from a run of ap_pipe.
"""


def main():
    script = sys.argv[0]
    try:
        repo = sys.argv[1]
    except IndexError:
        print('###')
        print('Run this script with 1-2 arguments: repo (required) and dbName (optional)')
        print('Please note that dbName path must be relative to repo')
        print('For example, python {0} rerun/MyAmazingRerun association.db'.format(script))
        print('###')
    else:
        rerunName = os.path.basename(os.path.normpath(repo))
        try:
            dbName = sys.argv[2]
        except IndexError:
            print('Using default dbName, association.db')
            print('Loading APDB...')
            objTable = loadAllApdbObjects(repo)
        else:
            print('Loading APDB...')
            objTable = loadAllApdbObjects(repo, dbName)
        finally:
            objFilter = setObjectFilter(objTable)
            print('Plotting Objects...')
            histPlot = plotDiaObjectHistogram(rerunName, objTable, objFilter)
            histName = rerunName + '_objHist.png'
            histPlot.savefig(histName)
            plt.show()
            if 'hits2015' in repo:
                skyPlot = plotDiaObjectsOnSky(rerunName, objTable, objFilter, hits=True)
            else:
                skyPlot = plotDiaObjectsOnSky(rerunName, objTable, objFilter, hits=False)
            skyName = rerunName + '_objSky.png'
            skyPlot.savefig(skyName)
            plt.show()
            print('Figures saved to {0} and {1}'.format(histName, skyName))


def connectToApdb(dbName, dbType='sqlite', schema=None):
    """Connect to an sqlite or postgres APDB.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, *full path* to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres

    Returns
    -------
    connection : `psycopg2.connect`
        A connection object to a database instance, ready for queries
    tables : `dict`
        Hard-wired tables known to exist in the APDB
    """
    if dbType == 'sqlite':
        connection = sqlite3.connect(dbName)
    elif dbType == 'postgres':
        if schema is None:
            raise RuntimeError('Schema must be set for postgres APDB')
        host = 'lsst-pg-devel1.ncsa.illinois.edu'
        connection = psycopg2.connect(dbname=dbName,
                                      host=host,
                                      options=f'-c search_path={schema}')
    else:
        raise ValueError('dbType must be sqlite or postgres')

    # These are some of the tables available in the APDB
    tables = {'obj': '"DiaObject"', 'src': '"DiaSource"',
              'sso': '"SSObject"', 'forcedsrc': '"DiaForcedSource"',
              'proto': '"ApdbProtoVisits"', 'match': '"DiaObject_To_Object_Match"'}

    return connection, tables


def loadAllApdbObjects(dbName, dbType='sqlite', schema=None):
    """Load a subset of DIAObject columns from a APDB into a pandas dataframe.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, *full path* to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres

    Returns
    -------
    objTable : `pandas.DataFrame`
        DIA Object Table containing only objects with validityEnd NULL.
        Columns selected are presently hard-wired here.
    """
    connection, tables = connectToApdb(dbName, dbType, schema)

    # Only get objects with validityEnd NULL because that means they are still valid
    objTable = pd.read_sql_query('select "diaObjectId", "ra", "decl", "nDiaSources", \
                                  "gPSFluxMean", "rPSFluxMean", "iPSFluxMean", \
                                  "zPSFluxMean", "yPSFluxMean", "validityEnd", "flags" from {0} \
                                  where "validityEnd" is NULL;'.format(tables['obj']), connection)
    return objTable

def loadAllApdbSources(dbName, dbType='sqlite', schema=None):
    """Load a subset of columns from all DIASources from a APDB into a pandas dataframe.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, full filepath to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Source Table including the columns hard-wired below.
    """
    connection, tables = connectToApdb(dbName, dbType, schema)

    # Load data from the source table
    srcTable = pd.read_sql_query('select "diaSourceId", "diaObjectId", \
                                  "ra", "decl", "ccdVisitId", \
                                  "midPointTai", "apFlux", "psFlux", "apFluxErr", \
                                  "psFluxErr", "totFlux", "totFluxErr", "x", "y", "snr", \
                                  "ixxPSF", "iyyPSF", "ixyPSF", "flags", "filterName" from {0}; \
                                  '.format(tables['src']), connection)
    return srcTable

def loadApdbSourcesByVisit(dbName, visit, dbType='sqlite', schema=None):
    """Load a subset of columns from all DIASources from a APDB into a pandas dataframe.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, full filepath to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres
    visit : `int`
        Visit number for loading objects

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Source Table including the columns hard-wired below.
    """
    connection, tables = connectToApdb(dbName, dbType, schema)
    # Load data from the source table
    srcTable = pd.read_sql_query('select "diaSourceId", "diaObjectId", \
                                  "ra", "decl", "ccdVisitId", \
                                  "midPointTai", "apFlux", "psFlux", "apFluxErr", \
                                  "psFluxErr", "totFlux", "totFluxErr", "x", "y", \
                                  "ixxPSF", "iyyPSF", "ixyPSF", "flags", "filterName" from {0} \
                                   where CAST("ccdVisitId" as text) like {1} ; \
                                  '.format(tables['src'], "'"+str(int(visit))+"%'"), connection)
    return srcTable

def loadApdbSourcesByBand(dbName, band, dbType='sqlite', schema=None):
    """Load a subset of columns from all DIASources from a APDB into a pandas dataframe.

    Parameters
    ----------
    dbName : `str`
        If dbType is sqlite, full filepath to the APDB on lsst-dev.
        If dbType is postgres, name of the APDB on lsst-pg-devel1.
    dbType : `str`, optional
        Either 'sqlite' or 'postgres'
    schema : `str`, optional
        Required if dbType is postgres
    band : `str`
        Band for loading objects (matched against filterName)

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Source Table including the columns hard-wired below.
    """
    connection, tables = connectToApdb(dbName, dbType, schema)
    # Load data from the source table
    srcTable = pd.read_sql_query('select "diaSourceId", "diaObjectId", \
                                  "ra", "decl", "ccdVisitId", \
                                  "midPointTai", "apFlux", "psFlux", "apFluxErr", \
                                  "psFluxErr", "totFlux", "totFluxErr", "x", "y", \
                                  "ixxPSF", "iyyPSF", "ixyPSF", "flags", "filterName" from {0} \
                                   where "filterName" = {1} ; \
                                  '.format(tables['src'], "'"+band+"'"), connection)
    return srcTable




def setObjectFilter(objTable):
    """Define a subset of objects to plot, i.e., make some kind of quality cut.

    The definition of objFilter is presently hard-wired here.

    Parameters
    ----------
    objTable : `pandas.DataFrame`
        DIA Object Table.

    Returns
    -------
    objFilter : `pandas.Series` of `bool`
        Filter applied to create a subset (e.g., quality cut) from objTable.
    """
    objFilter = ((objTable['nDiaSources'] > 14) & (objTable['flags'] == 0))
    numTotal = len(objTable['diaObjectId'])
    numFilter = len(objTable.loc[objFilter, 'diaObjectId'])
    print('There are {0} total DIAObjects and {1} filtered DIAObjects.'.format(numTotal, numFilter))
    return objFilter


def plotDiaObjectHistogram(objTable, objFiltered,
                           label1='All Objects', label2='Filtered Objects', title=''):
    """Create a histogram showing how many DIA Sources comprise the DIA Objects.

    Parameters
    ----------
    objTable : `pandas.DataFrame`
        DIA Object Table.
    objFiltered : `pandas.core.frame.DataFrame`
        DIA Object Table that is a filtered subset of objTable.
    label1 : `str`
        Legend label for the first DIA Object Table.
    label2 : `str`
        Legend label for the second (filtered) DIA Object Table.
    title : `str`
        Title for the plot, optional.

    Returns
    -------
    fig : `matplotlib.figure.Figure`
        Histogram of DIA Objects showing number of constituent DIA Sources.
    """
    fig = plt.figure()
    plt.xlabel('Number of Sources per Object', size=16)
    plt.ylabel('Object count', size=16)
    plt.ylim(0.7, 1e5)
    plt.yscale('log')
    binMax = np.max(objTable['nDiaSources'].values)
    plt.hist(objTable['nDiaSources'].values, bins=np.arange(0, binMax),
             color='#2979C1', label=label1)
    plt.hist(objFiltered['nDiaSources'].values, bins=np.arange(0, binMax),
             color='#Bee7F5', label=label2)
    plt.legend(frameon=False, fontsize=16)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title(title)
    return fig


def plotDiaObjectsOnSky(rerunName, objTable, objFilter, hits):
    """Create a plot of filtered DIAObjects on the sky.

    Parameters
    ----------
    rerunName : `str`
        Name of the directory at the end of the repo path.
    objTable : `pandas.DataFrame`
        DIA Object Table.
    objFilter : `pandas.Series` of `bool`
        Filter applied to create a subset (e.g., quality cut) from objTable.
    hits : `boolean`
        True for two panels with custom axes for the hits2015 dataset
        False for a single plot with automatic axis limits

    Returns
    -------
    fig : `matplotlib.figure.Figure`
        Scatter plot of the DIAObjects on the sky with RA, Dec axes.
        Point sizes and colors correspond to the number of sources per object.
    """
    fig = plt.figure(facecolor='white', figsize=(10, 8))

    if hits:  # two subplots
        dec_set1 = (objTable['decl'] > -2) & objFilter
        dec_set2 = (~dec_set1) & objFilter
        plt.subplots_adjust(wspace=0.1, hspace=0)

        # Panel 1: one HiTS field, on the right
        ax1 = plt.subplot2grid((100, 100), (0, 55), rowspan=90, colspan=45)
        plot1 = ax1.scatter(objTable.loc[dec_set1, 'ra'], objTable.loc[dec_set1, 'decl'],
                            marker='.', lw=0, s=objTable.loc[dec_set1, 'nDiaSources']*8,
                            c=objTable.loc[dec_set1, 'nDiaSources'], alpha=0.7,
                            cmap='viridis', linewidth=0.5, edgecolor='k')
        plt.xlabel('RA (deg)', size=16)
        plt.ylabel('Dec (deg)', size=16)
        ax1.spines['top'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position('right')
        ax1.invert_xaxis()  # RA should increase to the left
        loc = plticker.MultipleLocator(base=0.5)  # puts ticks at regular intervals
        ax1.yaxis.set_major_locator(loc)
        cb = fig.colorbar(plot1, orientation='horizontal', pad=0.3)
        cb.set_label('Number of Sources per Object', size=16)
        cb.set_clim(np.min(objTable.loc[dec_set2, 'nDiaSources']),
                    np.max(objTable.loc[dec_set2, 'nDiaSources']))
        cb.solids.set_edgecolor("face")
        cb.remove()  # don't show colorbar by this panel

        # Panel 2: two (overlapping) HiTS fields, on the left
        ax2 = plt.subplot2grid((100, 100), (0, 0), rowspan=90, colspan=50)
        plot2 = ax2.scatter(objTable.loc[dec_set2, 'ra'], objTable.loc[dec_set2, 'decl'],
                            marker='.', lw=0, s=objTable.loc[dec_set2, 'nDiaSources']*8,
                            c=objTable.loc[dec_set2, 'nDiaSources'], alpha=0.7,
                            cmap='viridis', linewidth=0.5, edgecolor='k')
        plt.xlabel('RA (deg)', size=16)
        plt.ylabel('Dec (deg)', size=16)
        plt.title(rerunName, size=16)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.invert_xaxis()
        cax = plt.subplot2grid((100, 100), (60, 55), rowspan=5, colspan=45)
        cb2 = fig.colorbar(plot2, cax=cax, orientation='horizontal', pad=0.1)
        cb2.set_label('Number of Sources per Object', size=16)
        cb2.set_clim(np.min(objTable.loc[dec_set2, 'nDiaSources']),
                     np.max(objTable.loc[dec_set2, 'nDiaSources']))
        cb2.solids.set_edgecolor("face")

    else:  # one main plot
        ax = fig.add_subplot(111)
        plot = ax.scatter(objTable.loc[objFilter, 'ra'],
                          objTable.loc[objFilter, 'decl'], marker='.', lw=0,
                          s=objTable.loc[objFilter, 'nDiaSources']*8,
                          c=objTable.loc[objFilter, 'nDiaSources'],
                          alpha=0.7, cmap='viridis', linewidth=0.5, edgecolor='k')
        plt.xlabel('RA (deg)', size=16)
        plt.ylabel('Dec (deg)', size=16)
        plt.title(rerunName, size=16)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.invert_xaxis()  # RA should increase to the left
        cb = fig.colorbar(plot, orientation='horizontal')
        cb.set_label('Number of Sources per Object', size=16)
        cb.set_clim(np.min(objTable.loc[objFilter, 'nDiaSources']),
                    np.max(objTable.loc[objFilter, 'nDiaSources']))
        cb.solids.set_edgecolor("face")
    return fig


if __name__ == '__main__':
    main()
