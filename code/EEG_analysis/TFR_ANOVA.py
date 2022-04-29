#   This script is part (5) of the EEG analysis of an antisaccade task and 
#   creates time-frequency representations of the preprocessed & epoched data,
#   then runs ANOVA to compared conditions and creates figures 

#  Copyright (C) April 2020, last modified April 2022
#   J.Waldthaler, A. Sperlich, D. Pedrosa
#   University Hospital of Gießen and Marburg
#
#   This software may be used, copied, or redistributed as long as it is
#   not sold and this copyright notice is reproduced on each copy made.
#   This routine is provided as is without any express or implied
#   warranties whatsoever.


import os.path as op
import numpy as np
import mne
from mne.epochs import equalize_epoch_counts
import matplotlib.pyplot as plt
from mne.stats import (f_threshold_mway_rm,f_mway_rm, permutation_cluster_1samp_test)
from scipy import stats
from config import subjects, eeg_path, sub_path

# define ROI
# list of channels for lat prefrontal ROI:
channels = {
    'F6',
    'F8',
    'AF8',
    'FC6'
}
# list of channels for midfrontal ROI:
channels = {
    'Fz',
    'F1',
    'F2'
}

freqs = np.arange(4, 30, 1)
cycles = freqs/2
sfreq=250
n_freqs= 26
tmin, tmax = 1.0, 2.0 # 1.0 to 2.0 = complete period of fixation cross presentation, epochs range from 0.5 to 2.1
n_times=int((((tmax-tmin)*1000))/4+1)
n_subjects=len(subjects)
baseline = (0.8, 0.9) # ISI part of epoch = 0.5 to 1.0
p_threshold = 0.05
n_permutations = 1000

#for ANOVA:
n_conditions =  3
factor_levels = [3]  # number of levels in each factor
effects = 'A'
f_threshold = f_threshold_mway_rm(n_subjects, factor_levels, effects, p_threshold)
tail = 1  # f-test, so tail > 0

# for pairwise permutation t tests
t_threshold = stats.distributions.t.ppf(1 - p_threshold / 2, n_subjects - 1)

# create arrays to append and store data of all subjects
X_60 = np.zeros((0, n_freqs, n_times))
X_off = np.zeros((0, n_freqs, n_times))
X_130 = np.zeros((0, n_freqs, n_times))

# loop TFR for all subjects
for subj in subjects:
    print('Now processing subject ' + subj)
    sub_path = op.join(eeg_path, subj)
    epochs_off_fname = op.join(sub_path, subj + '_off_epo.fif')
    epochs_130_fname = op.join(sub_path, subj + '_130_epo.fif')
    epochs_60_fname = op.join(sub_path, subj + '_60_epo.fif')

    # read data
    epo_off = mne.read_epochs(epochs_off_fname)
    epo_off.pick_channels(ch_names=channels)
    epo_130 = mne.read_epochs(epochs_130_fname)
    epo_130.pick_channels(ch_names=channels)
    epo_60 = mne.read_epochs(epochs_60_fname)
    epo_60.pick_channels(ch_names=channels)

    # select only correct trials
    epo_off =   epo_off["1"]
    epo_130 =   epo_130["1"]
    epo_60  =   epo_60["1"]

    # equalize nr. of epochs for all three conditions for each subject, truncate = take epochs from the ends of the files
    equalize_epoch_counts([epo_off, epo_130, epo_60], method="truncate")
    # calculate TFR with baseline correction and cropping for each condition
    tfr_off = mne.time_frequency.tfr_morlet(epo_off, freqs=freqs, n_cycles=cycles, n_jobs=1, average=True,
                                           use_fft=True, return_itc=False)
    tfr_off = tfr_off.apply_baseline(baseline, 'logratio')
    tfr_off.crop(tmin, tmax)

    tfr_130 = mne.time_frequency.tfr_morlet(epo_130, freqs=freqs, n_cycles=cycles, n_jobs=1, average=True,
                                           use_fft=True, return_itc=False)
    tfr_130 = tfr_130.apply_baseline(baseline, 'logratio')
    tfr_130.crop(tmin, tmax)

    tfr_60 = mne.time_frequency.tfr_morlet(epo_60, freqs=freqs, n_cycles=cycles, n_jobs=1, average=True,
                                           use_fft=True, return_itc=False)
    tfr_60=tfr_60.apply_baseline(baseline, 'logratio')
    tfr_60.crop(tmin, tmax)

    # get np arrays of epochs data
    epo_off = tfr_off.data
    epo_130 = tfr_130.data
    epo_60 = tfr_60.data

    #average across channels to create one observation per ROI
    epo_off = np.mean(epo_off, axis=0)
    epo_130 = np.mean(epo_130, axis=0)
    epo_60 = np.mean(epo_60, axis=0)
    epo_off = np.expand_dims(epo_off, 0)
    epo_130 = np.expand_dims(epo_130, 0)
    epo_60 = np.expand_dims(epo_60, 0)

    # append this subject's data to array
    X_off = np.append(X_off, epo_off, axis=0)
    X_130 = np.append(X_130, epo_130, axis=0)
    X_60 = np.append(X_60, epo_60, axis=0)

# average across subjects for figures
X_off_mean = np.mean(X_off, axis=0)
X_130_mean = np.mean(X_130, axis=0)
X_60_mean = np.mean(X_60, axis=0)

# plot all conditions in one figure
vmax = 0.2
cmap='RdYlBu_r'
cmap2='Reds'
times = 1e3 * tfr_60.times -2000 # change unit to ms and adjust to 0 ms = target onset
fig, (plt1, plt2, plt3, cax) = plt.subplots(nrows=1, ncols=4, figsize=(12, 3), gridspec_kw={"width_ratios":[1,1,1, 0.05]})
fig.subplots_adjust(wspace=0.15)
plt1a = plt1.imshow((X_off_mean),extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap,
           vmin=-vmax, vmax=vmax)
plt1.set_xlabel('Time (ms)')
plt1.set_ylabel('Frequency (Hz)')
plt1.set_title('DBS off')
plt2.imshow(X_130_mean,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap,
           vmax=vmax, vmin=-vmax)
plt2.set_xlabel('Time (ms)')
plt2.set_title('130 Hz DBS')
plt2.set_yticklabels([])
plt3.imshow(X_60_mean,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap,
           vmax=vmax, vmin=-vmax)
plt3.set_yticklabels([])
plt3.set_xlabel('Time (ms)')
plt3.set_title('60 Hz DBS')
fig.colorbar(plt1a, cax=cax)
cax.set_ylabel('Power change from baseline (dB)')
plt.tight_layout()

###########################STATS#################################
# RM-ANOVA
# create new tuple to store all three conditions, each as single array
X_1 = (X_off, X_130, X_60)

# define the ANOVA statistic functions
# needs swaping of axes since ANOVA expects an input array of dimensions: subjects X conditions X observations
def stat_fun(*args):
     return f_mway_rm(np.swapaxes(args, 1,0), factor_levels=factor_levels, effects=effects, return_pvals=False)[0]
# run ANOVA permutation test with 1000 permutations
T_obs, clusters, cluster_p_values, h0 = mne.stats.permutation_cluster_test(X_1, stat_fun=stat_fun, threshold=f_threshold,
                                                                            n_jobs=1, tail=tail,
                                                                           n_permutations=n_permutations, buffer_size=None,
                                                                           out_type='mask')
# check at p values of all clusters
print(cluster_p_values)
# select only significant clusters for figure
good_clusters = np.where(cluster_p_values < p_threshold)[0]
T_obs_plot = T_obs.copy()
T_obs_plot[~clusters[np.squeeze(good_clusters)]] = np.nan

# PAIRWISE COMPARISONS WITH PERMUTATION T TESTS
# create new array with data from all subjects and all conditions
X = np.zeros((n_subjects, n_freqs, n_times, 3))
X[:, :, :,0] += X_off
X[:, :, :,1] += X_130
X[:, :, :,2] += X_60
# make paired contrast
Con_1_0 = X[:, :, :, 1] - X[:, :, :, 0]
Con_2_0 = X[:, :, :, 2] - X[:, :, :, 0]
Con_2_1 = X[:, :, :, 2] - X[:, :, :, 1]
# average across subjects
Con_1_0_mean = np.mean(Con_1_0, axis=0)
Con_2_0_mean = np.mean(Con_2_0, axis=0)
Con_2_1_mean = np.mean(Con_2_1, axis=0)

# calculate permutation tests
# 130 Hz vs off
T_obs_1_0, clusters_1_0, cluster_p_values_1_0, H0 = \
    permutation_cluster_1samp_test(Con_1_0, n_permutations=n_permutations, threshold=t_threshold)
# see p values of clusters
print(cluster_p_values_1_0)
# Create stats image with only significant clusters
T_obs_plot_1_0 = np.nan  * np.ones_like(T_obs_1_0)
for c, p_val in zip(clusters_1_0, cluster_p_values_1_0):
   if p_val <= p_threshold:
     T_obs_plot_1_0[c] = T_obs_1_0[c]

# 60 Hz vs. off
T_obs_2_0, clusters_2_0, cluster_p_values_2_0, H0 = \
    permutation_cluster_1samp_test(Con_2_0, n_permutations=n_permutations)
# see p values of clusters
print(cluster_p_values_2_0)
# Create stats image with only significant clusters
T_obs_plot_2_0 = np.nan  * np.ones_like(T_obs_2_0)
for c, p_val in zip(clusters_2_0, cluster_p_values_2_0):
   if p_val <= p_threshold:
     T_obs_plot_2_0[c] = T_obs_2_0[c]

# 60 Hz vs. 130 Hz
T_obs_2_1, clusters_2_1, cluster_p_values_2_1, H0 = \
    permutation_cluster_1samp_test(Con_2_1, n_permutations=n_permutations)
# see p values of clusters
print(cluster_p_values_2_1)
# Create stats image with only significant clusters
T_obs_plot_2_1 = np.nan  * np.ones_like(T_obs_2_1)
for c, p_val in zip(clusters_2_1, cluster_p_values_2_1):
   if p_val <= p_threshold:
     T_obs_plot_2_1[c] = T_obs_2_1[c]

# plot ANOVA and pairwise contrasts of conditions in one figure
vmax = 0.25
cmap3 = "RdBu_r"
cmap2="gist_heat_r"

vmin2 = -3
vmax2= 3
fig, (plt4, cax2, cax1) = plt.subplots(1,3, figsize=(6, 4), gridspec_kw={"width_ratios":[1, 0.05, 0.05]})

plt4b=plt4.imshow(T_obs,  extent=[times[0], times[-1],freqs[0], freqs[-1]], aspect='auto',origin='lower', cmap='binary_r')
plt4a=plt4.imshow(T_obs_plot,  extent=[times[0], times[-1],freqs[0], freqs[-1]], aspect='auto',origin='lower', cmap=cmap2)
plt4.set_xlabel('Time (ms)')
plt4.set_ylabel('Frequency (Hz)')
plt4.set_title('RM-ANOVA')
fig.colorbar(plt4a, cax=cax1)
cax1.set_ylabel('F values of significant clusters\n')
fig.colorbar(plt4b, cax=cax2)
cax2.set_ylabel('F values\n')
plt.tight_layout(w_pad=0.1)

fig, (plt1, plt2,  plt3, cax3) = plt.subplots(1,4, figsize=(15, 4),
                                                               gridspec_kw={"width_ratios":[1,1, 1,0.05]})
plt1a = plt1.imshow(Con_1_0_mean,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower',
                    cmap=cmap,vmin=-vmax, vmax=vmax, alpha=0.75)
plt1b=plt1.imshow(T_obs_plot_1_0,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap, vmax=vmax2, vmin = vmin2)
plt1.set_title('130 Hz - DBS off Contrast')
plt1.set_xlabel('Time (ms)')
plt1.set_ylabel('Frequency (Hz)')
plt2.imshow(Con_2_0_mean,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap,
           vmax=vmax, vmin=-vmax, alpha = 0.75)
plt2b = plt2.imshow(T_obs_plot_2_0,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap, vmax = vmax2, vmin= vmin2)
plt2.set_title('60 Hz - DBS off Contrast')
plt2.set_yticklabels([])
plt3.imshow(Con_2_1_mean,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap,
           vmax=vmax, vmin=-vmax, alpha = 0.75)
plt3.imshow(T_obs_plot_2_1,extent=[times[0], times[-1], freqs[0], freqs[-1]],aspect='auto', origin='lower', cmap=cmap)
plt3.set_yticklabels([])
plt3.set_title('60 Hz - 130 Hz Contrast')
fig.colorbar(plt1a,cax=cax3)
cax3.set_ylabel('Difference in power \n change from baseline (dB)\n ')
#plt.tight_layout(w_pad=0.1)

X_off_meanmean = np.mean(X_off_mean, axis=0)
X_130_meanmean = np.mean(X_130_mean, axis=0)
X_60_meanmean = np.mean(X_60_mean, axis=0)
plt.figure(figsize=(6, 4))
plt.plot(times, X_off_meanmean, label='off-DBS', color='black')
plt.plot(times, X_130_meanmean, label='130 Hz DBS', color='red')
plt.plot(times, X_60_meanmean, label='60 Hz DBS', color='blue')
plt.ylim(-0.10, 0.10)
plt.xlabel('Time (ms)')
plt.ylabel('Power change from baseline (dB)')
plt.legend()
plt.title('Beta power (18-26 Hz)', fontsize=12)
plt.show()
