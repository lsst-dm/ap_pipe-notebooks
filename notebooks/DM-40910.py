import lsst.daf.butler as dafButler
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from mask_tools import should_ignore, get_enabled_flags_names
import numpy as np

butler = dafButler.Butler("/dc2/dc2")

# diff_collection = 'u/nsedaghat/DLDataset0.1.1-diffs/20220522T064105Z'
# calexp_collection = 'u/nsedaghat/Dataset0.1.1-sciCalexps/20220521T042355Z'

# diff_collection = 'u/nsedaghat/DLDataset0.1.2-diffs/20231004T201903Z'
# calexp_collection = 'u/nsedaghat/RegeneratedCalexps_w_2022_20_nobfk/20231002T214649Z'

diff_collection = 'u/nsedaghat/DLDataset0.1.2-diffs_w_2023_38/20231007T024858Z'
calexp_collection = 'u/nsedaghat/RegeneratedCalexps_w_2022_20_nobfk/20231002T214649Z'

def xmatch_det_truth(df_true, df_det,
                     ra_col_true="ra", dec_col_true="dec",
                     ra_col_det="ra", dec_col_det="decl"
                     ):
    if len(df_det) == 0:
        return []
    elif len(df_true) == 0:
        return False * len(df_det)

    cat1 = SkyCoord(ra=df_true[ra_col_true].values * u.degree,
                    dec=df_true[dec_col_true].values * u.degree)

    cat2 = SkyCoord(ra=df_det[ra_col_det].values * u.degree,
                    dec=df_det[dec_col_det].values * u.degree)

    idx, d2d, d3d = cat2.match_to_catalog_sky(cat1)

    imatch = d2d < 1 * u.arcsec  # .0005*u.degree
    return imatch


def run_single(i):
    ''' Extracts the truth and detection catalogs for a single visit.

    Parameters
    ----------
    i : int
        The index of the visit to extract.
        This the index in the list of visits returned by the butler query,
        not the visit ID.
    '''


    dataId = ids[i]
    print(dataId)

    # Load the truth (by fetching the calexp's URI)
    sciFile = butler.getURI('calexp',
                            collections=[calexp_collection],
                            dataId=dataId).geturl().replace('file://', '')

    truthfile = sciFile + '.SNVAR.parq'
    truth = pd.read_parquet(truthfile, engine='pyarrow')

    # Load the source table
    src_table = butler.get('goodSeeingDiff_diaSrcTable',
                           dataId=dataId,
                           collections=[diff_collection])


    #--- Obtain mask info for truth rows, and apply
    sci = butler.get('calexp',
                     collections=[calexp_collection],
                     dataId=dataId)
    X,Y = sci.getWcs().skyToPixelArray(ra=truth['ra'],dec=truth['dec'],degrees=True)
    truth = truth.assign(x=X,y=Y) # add x,y columns to truth
    mask = sci.getMaskedImage().getMask()
    truth['mask'] = mask.getArray()[truth['y'].astype(int),truth['x'].astype(int)]

    truth = truth[~truth.apply(lambda row: should_ignore(row['mask']), axis=1)]


    xmatch_result = xmatch_det_truth(truth, src_table, dec_col_det='dec')
    src_table['label'] = xmatch_result

    return src_table, truth

ids = list(butler.registry.queryDataIds(['visit','detector'],
                                   datasets='goodSeeingDiff_diaSrcTable',
                                        collections=[diff_collection],
                                        #where="instrument='LSSTCam-imSim' AND skymap='DC2' AND exposure<1100000"
                                        ))


all_src = pd.DataFrame()
all_truth = pd.DataFrame()
for i in range(min(20, len(ids))):
    src_table, truth = run_single(i)
    all_src = pd.concat([all_src, src_table])
    all_truth = pd.concat([all_truth, truth])

# The "label" column is True for all sources that are matched to a truth object.
# Let's look at some stats.
print(f"Found {len(all_src[all_src['label']])} matches out of {len(all_src)} sources")

## Now sweep over the SNR threshold and create a ROC curve
import matplotlib.pyplot as plt
max_snr = all_src['snr'].max()
snr_thresholds = np.linspace(0, max_snr, 100)
true_positives = []
false_positives = []
for snr_threshold in snr_thresholds:
    true_positives.append(len(all_src[(all_src['label']) & (all_src['snr'] > snr_threshold)]))
    false_positives.append(len(all_src[(~all_src['label']) & (all_src['snr'] > snr_threshold)]))

plt.plot(false_positives, true_positives, label='ROC curve')
plt.plot([0, len(all_src)], [0, len(all_src)], linestyle='--', label='Random')
plt.xlabel("False positives")
plt.ylabel("True positives")
plt.title("ROC curve")
plt.legend()
plt.show()

## and a Precision-Recall curve
precision = np.array(true_positives) / (np.array(true_positives) + np.array(false_positives) + 1e-10)
recall = np.array(true_positives) / len(all_truth)
plt.plot(recall, precision, label='Precision-Recall curve')
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall curve")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.show()
