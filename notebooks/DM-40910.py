##
import lsst.daf.butler as dafButler
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from mask_tools import should_ignore, get_enabled_flags_names
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

butler = dafButler.Butler("/dc2/dc2")

# diff_collection = 'u/nsedaghat/DLDataset0.1.1-diffs/20220522T064105Z'
# calexp_collection = 'u/nsedaghat/Dataset0.1.1-sciCalexps/20220521T042355Z'

# diff_collection = 'u/nsedaghat/DLDataset0.1.2-diffs/20231004T201903Z'
# calexp_collection = 'u/nsedaghat/RegeneratedCalexps_w_2022_20_nobfk/20231002T214649Z'

diff_collection = 'u/nsedaghat/DLDataset0.1.2-diffs_w_2023_38/20231007T024858Z'
calexp_collection = 'u/nsedaghat/RegeneratedCalexps_w_2022_20_nobfk/20231002T214649Z'

def xmatch_det_truth(df_true, df_det,
                     ra_col_true="ra", dec_col_true="dec",
                     ra_col_det="ra", dec_col_det="dec"
                     ):
    if len(df_det) == 0:
        return []
    elif len(df_true) == 0:
        return [-1] * len(df_det)

    cat1 = SkyCoord(ra=df_true[ra_col_true].values * u.degree,
                    dec=df_true[dec_col_true].values * u.degree)

    cat2 = SkyCoord(ra=df_det[ra_col_det].values * u.degree,
                    dec=df_det[dec_col_det].values * u.degree)

    idx, d2d, d3d = cat2.match_to_catalog_sky(cat1)

    imatch_det = d2d < 1 * u.arcsec  # .0005*u.degree

   # For matched detections, store the truth index. Otherwise, store -1.
    matched_indices = np.where(d2d < 1 * u.arcsec, idx, -1)

    return matched_indices.tolist()

# # Create a boolean array for the truth table
#     imatch_true = np.zeros(len(df_true), dtype=bool)
#     imatch_true[idx[imatch_det]] = True

#     return imatch_det.tolist(), imatch_true.tolist()

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

    # Cross-match the truth and detections
    imatch_det = xmatch_det_truth(truth, src_table)
    src_table['matched_truth_index'] = imatch_det
    mask = src_table['matched_truth_index'] != -1
    src_table.loc[mask, 'matched_truth_magr'] = truth['mag_r'].values[src_table.loc[mask, 'matched_truth_index'].values]
    #src_table['matched_truth_magr'] = truth['mag_r'].values[src_table['matched_truth_index'].values]
    src_table['label'] = src_table['matched_truth_index'] != -1
    truth['detected'] = truth.index.isin(imatch_det)

    # imatch_det, imatch_true = xmatch_det_truth(truth, src_table, dec_col_det='dec')
    # src_table['label'] = imatch_det
    # truth['detected'] = imatch_true

    return src_table, truth

ids = list(butler.registry.queryDataIds(['visit','detector'],
                                   datasets='goodSeeingDiff_diaSrcTable',
                                        collections=[diff_collection],
                                        #where="instrument='LSSTCam-imSim' AND skymap='DC2' AND exposure<1100000"
                                        ))

## print the data ids in a readable format
N_test = 100
for dataId in ids[:N_test]:
    print(dataId.to_simple())

## Loop on data ids and extract the truth and detection catalogs
all_src = pd.DataFrame()
all_truth = pd.DataFrame()

for i in range(min(len(ids), N_test)):
    src_table, truth = run_single(i)
    all_src = pd.concat([all_src, src_table])
    all_truth = pd.concat([all_truth, truth])

# The "label" column is True for all sources that are matched to a truth object.
# Let's look at some stats.
print("Total number of detected sources: ", len(all_src))
print("Number of sources matching the truth: ", np.sum(all_src['label']))
print("Number of sources not matching the truth (false positives): ", np.sum(1-all_src['label']))
print("Total number of truth objects: ", len(all_truth))
print("Number of truth objects not matched to a source (misses): ", len(all_truth) - np.sum(all_src['label']))

## Some statistics about detections
# plot the hist of mag_r of true objects that are detected.
# you need to drop the NaNs and infinities first.
plt.figure(figsize=(10, 8))
plt.hist(all_truth['mag_r'][np.isfinite(all_truth['mag_r']) & (all_truth['detected'] == False)], bins=100, label='missed')
plt.hist(all_truth['mag_r'][np.isfinite(all_truth['mag_r']) & all_truth['detected']], bins=100, alpha=0.5, label='detected')
plt.xlabel('mag_r')
plt.ylabel('Number of true objects')
plt.legend()
plt.savefig('mag_r_hist.png')
plt.show()


## make the same histograms plot, but with delta_flux
plt.figure(figsize=(10, 8))
plt.hist(all_truth['delta_flux'][np.isfinite(all_truth['delta_flux']) & (all_truth['detected'] == False)],
         bins=100, range=(-5000, 5000), label='missed')
plt.hist(all_truth['delta_flux'][np.isfinite(all_truth['delta_flux']) & all_truth['detected']],
         bins=100, range=(-5000, 5000), alpha=0.5, label='detected')
plt.xlabel('delta_flux')
plt.ylabel('Number of true objects')
plt.legend()
plt.savefig('delta_flux_hist.png')
plt.show()


## Compute the true and false posititves.
def eval_filtered(truth_magr_max=np.inf):
    ''' Computes the true and false positives for a given maximum truth mag_r
    '''

    good_sources = all_src[all_src['label'] & (all_src['matched_truth_magr'] < truth_magr_max)]

    max_snr = all_src['snr'].max()
    snr_thresholds = np.linspace(0, max_snr, 100)
    true_positives = []
    false_positives = []
    for snr_threshold in snr_thresholds:
        true_positives.append(np.sum(good_sources['snr'] > snr_threshold))
        false_positives.append(len(all_src[(~all_src['label']) & (all_src['snr'] > snr_threshold)]))

    true_positives = np.array(true_positives)
    false_positives = np.array(false_positives)
    return true_positives, false_positives, snr_thresholds

## The ROC curve
plt.figure(figsize=(10, 10))
for mag_r in range(10, 30, 2):
    true_positives, false_positives, snr_thresholds = eval_filtered(truth_magr_max=mag_r)
    plt.plot(false_positives/len(all_src), true_positives/len(all_src), label='mag_r < {}'.format(mag_r))

plt.plot([0, 1], [0, 1], linestyle='--', label='Reference')
#plt.plot([0, len(all_src)], [0, len(all_src)], linestyle='--', label='Reference')
#plt.xlim(0, len(all_src))
plt.xlabel("False positives")
plt.ylabel("True positives")
plt.title("ROC curve")
plt.legend()
plt.savefig("roc_curve.png")
plt.show()

## ROC curve using scikit-learn
from sklearn.metrics import roc_curve, auc

fpr, tpr, thresholds = roc_curve(all_src['label'], all_src['snr'])

## and a set of Precision-Recall curve
plt.figure(figsize=(10, 10))
for mag_r in range(10, 30, 2):
    true_positives, false_positives, snr_thresholds = eval_filtered(mag_r)
    precision = np.array(true_positives) / (np.array(true_positives) + np.array(false_positives) + 1e-10)
    recall = np.array(true_positives) / len(all_truth)
    plt.plot(recall, precision, label=f'mag_r < {mag_r}', linewidth=5, alpha=0.5)

plt.xlabel("Recall/Completeness")
plt.ylabel("Precision/Purity")
plt.title("Precision-Recall curve")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.savefig("precision_recall_curve.png")
plt.show()

## And a curve showing completeness and purity vs. SNR threshold
plt.figure(figsize=(20, 10))
plt.subplot(1, 2, 1)
for mag_r in range(10, 30, 2):
    true_positives, false_positives, snr_thresholds = eval_filtered(mag_r)
    plt.plot(snr_thresholds, true_positives/len(all_truth), label=f'mag_r < {mag_r}', linewidth=5, alpha=0.5)

plt.xlim(0, max(snr_thresholds))
plt.ylim(0, 1)
plt.xlabel("SNR threshold")
# use latex
plt.ylabel("Completeness: $TP/(TP+FN)$")
plt.legend()

plt.subplot(1, 2, 2)
for mag_r in range(10, 30, 2):
    true_positives, false_positives, snr_thresholds = eval_filtered(mag_r)
    plt.plot(snr_thresholds, true_positives/(true_positives + false_positives + 1e-10), label=f'mag_r < {mag_r}', linewidth=5, alpha=0.5)

plt.xlim(0, max(snr_thresholds))
plt.ylim(0, 1)
plt.xlabel("SNR threshold")
plt.ylabel("Purity: $TP/(TP+FP)$")
plt.legend()


plt.savefig("detection_rate.png")
plt.show()
