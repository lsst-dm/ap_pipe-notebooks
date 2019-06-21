import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
import sqlite3

"""Script to make plots of DIAObjects using a Prompt Products Database
(PPDB) resulting from a run of ap_pipe.
"""


def main():
    script = sys.argv[0]
    try:
        repo = sys.argv[1]
    except:
        print('###')
        print('Run this script with 1-2 arguments: repo (required) and dbName (optional)')
        print('Please note that dbName path must be relative to repo')
        print('For example, python {0} rerun/MyAmazingRerun association.db'.format(script))
        print('###')
    else:
        rerunName = os.path.basename(os.path.normpath(repo))
        try:
            dbName = sys.argv[2]
        except:
            print('Using default dbName, association.db')
            print('Loading PPDB...')
            objTable = loadAllPpdbObjects(repo)
        else:
            print('Loading PPDB...')
            objTable = loadAllPpdbObjects(repo, dbName)
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


def loadAllPpdbObjects(repo, dbName='association.db'):
    """Load select DIAObject columns from a PPDB into a pandas dataframe.

    Parameters
    ----------
    repo : `str`
        Path to an output repository from an ap_pipe run.
    dbName : `str`, optional
        Name of the PPDB, which must reside in (or relative to) repo.

    Returns
    -------
    objTable : `pandas.DataFrame`
        DIA Object Table containing only objects with validityEnd NULL.
        Columns selected are presently hard-wired here.
    """
    connection = sqlite3.connect(os.path.join(repo, dbName))

    # These are some of the tables available in the ppdb
    tables = {'obj': 'DiaObject', 'src': 'DiaSource'}

    # Only get objects with validityEnd NULL because that means they are still valid
    objTable = pd.read_sql_query('select diaObjectId, ra, decl, nDiaSources, \
                                  gPSFluxMean, validityEnd, flags from {0} \
                                  where validityEnd is NULL;'.format(tables['obj']), connection)
    return objTable


def loadAllPpdbSources(repo, dbName='association.db'):
    """Load select columns from all DIASources from a PPDB into a pandas dataframe.

    Parameters
    ----------
    repo : `str`
        Path to an output repository from an ap_pipe run.
    dbName : `str`, optional
        Name of the PPDB, which must reside in (or relative to) repo.
        dbPath):

    Returns
    -------
    srcTable : `pandas.DataFrame`
        DIA Source Table including the columns hard-wired below.
    """
    connection = sqlite3.connect(os.path.join(repo, dbName))

    # These are some of the tables available in the ppdb
    tables = {'obj': 'DiaObject', 'src': 'DiaSource'}

    # Load data from the source table
    srcTable = pd.read_sql_query('select diaSourceId, diaObjectId, \
                                  ra, decl, ccdVisitId, \
                                  midPointTai, apFlux, psFlux, apFluxErr, \
                                  psFluxErr, totFlux, totFluxErr, flags from {0} \
                                  '.format(tables['src']), connection)
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


def plotDiaObjectHistogram(rerunName, objTable, objFilter):
    """Create a histogram showing all DIAObjects and filtered DIAObjects.

    Parameters
    ----------
    rerunName : `str`
        Name of the directory at the end of the repo path.
    objTable : `pandas.DataFrame`
        DIA Object Table.
    objFilter : `pandas.Series` of `bool`
        Filter applied to create a subset (e.g., quality cut) from objTable.

    Returns
    -------
    fig : `matplotlib.figure.Figure`
        Histogram of all DIAObjects in objTable with a histogram of the
        subset described by objFilter plotted on top.
    """
    fig = plt.figure()
    plt.xlabel('Number of Sources per Object', size=16)
    plt.ylabel('Object count', size=16)
    plt.title(rerunName, size=16)
    plt.ylim(0.7, 1e5)
    plt.yscale('log')
    binMax = np.max(objTable['nDiaSources'].values)
    plt.hist(objTable['nDiaSources'].values, bins=np.arange(0, binMax),
             color='#2979C1', label='All Objects')
    plt.hist(objTable.loc[objFilter, 'nDiaSources'].values, bins=np.arange(0, binMax),
             color='#Bee7F5', label='Filtered Objects')
    plt.legend(frameon=False, fontsize=16)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
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
