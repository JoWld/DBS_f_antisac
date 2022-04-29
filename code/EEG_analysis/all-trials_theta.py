#   This script is part of the EEG analysis of an antisaccade task
#   the output is an Excel file with theta power changes in dB for all trials

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
import pandas as pd
from config import subjects, eeg_path, sub_path


file_name = 'all-trials_theta_reg.xlsx'

# list of channels for midfrontal ROI:
channels = {
    'Fz',
    'F1',
    'F2'
}

freqs = np.arange(3, 31, 1)
cycles = freqs/2
sfreq=250

# select frequencies, 4-8 for theta
fmin = 4
fmax = 8

tmin, tmax = 1.6, 2.0 
baseline = (0.8, 0.9)

# create a dataframe to store data
all = pd.DataFrame(columns=['subj', 'theta'])

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

    # select correct trials
    epo_off_c =   epo_off["1"]
    epo_130_c =   epo_130["1"]
    epo_60_c  =   epo_60["1"]

    # select error trials
    epo_off_e =   epo_off["2"]
    epo_130_e =   epo_130["2"]
    epo_60_e  =   epo_60["2"]

    # calculate TFR with baseline correction and cropping for each condition
    tfr_off_c = mne.time_frequency.tfr_morlet(epo_off_c, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_off_c.apply_baseline(baseline, 'logratio')
    tfr_off_c.crop(tmin, tmax, fmin, fmax)

    tfr_130_c = mne.time_frequency.tfr_morlet(epo_130_c, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_130_c.apply_baseline(baseline, 'logratio')
    tfr_130_c.crop(tmin, tmax, fmin, fmax)

    tfr_60_c = mne.time_frequency.tfr_morlet(epo_60_c, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_60_c.apply_baseline(baseline, 'logratio')
    tfr_60_c.crop(tmin, tmax, fmin, fmax)

    tfr_off_e = mne.time_frequency.tfr_morlet(epo_off_e, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_off_e.apply_baseline(baseline, 'logratio')
    tfr_off_e.crop(tmin, tmax, fmin, fmax)

    tfr_130_e = mne.time_frequency.tfr_morlet(epo_130_e, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_130_e.apply_baseline(baseline, 'logratio')
    tfr_130_e.crop(tmin, tmax, fmin, fmax)

    tfr_60_e = mne.time_frequency.tfr_morlet(epo_60_e, freqs=freqs, n_cycles=cycles, n_jobs=1, average=False,
                                           use_fft=True, return_itc=False)
    tfr_60_e.apply_baseline(baseline, 'logratio')
    tfr_60_e.crop(tmin, tmax, fmin, fmax)

    # get np arrays of epochs data
    epo_off_c = tfr_off_c.data[:, :, :]
    epo_130_c = tfr_130_c.data[:, :, :]
    epo_60_c = tfr_60_c.data[:, :, :]
    epo_off_e = tfr_off_e.data[:, :, :]
    epo_130_e = tfr_130_e.data[:, :, :]
    epo_60_e = tfr_60_e.data[:, :, :]

    #average across channels and frequencies to create one observation per ROI
    epo_off_c = np.mean(epo_off_c, axis=1)
    epo_130_c = np.mean(epo_130_c, axis=1)
    epo_60_c = np.mean(epo_60_c, axis=1)
    epo_off_c = np.mean(epo_off_c, axis=1)
    epo_130_c = np.mean(epo_130_c, axis=1)
    epo_60_c = np.mean(epo_60_c, axis=1)
    epo_off_c = np.mean(epo_off_c, axis=1)
    epo_130_c = np.mean(epo_130_c, axis=1)
    epo_60_c = np.mean(epo_60_c, axis=1)
    epo_off_e = np.mean(epo_off_e, axis=1)
    epo_130_e = np.mean(epo_130_e, axis=1)
    epo_60_e = np.mean(epo_60_e, axis=1)
    epo_off_e = np.mean(epo_off_e, axis=1)
    epo_130_e = np.mean(epo_130_e, axis=1)
    epo_60_e = np.mean(epo_60_e, axis=1)
    epo_off_e = np.mean(epo_off_e, axis=1)
    epo_130_e = np.mean(epo_130_e, axis=1)
    epo_60_e = np.mean(epo_60_e, axis=1)

    epo_60_e = epo_60_e.tolist()
    epo_60_c = epo_60_c.tolist()
    epo_130_e = epo_130_e.tolist()
    epo_130_c = epo_130_c.tolist()
    epo_off_e = epo_off_e.tolist()
    epo_off_c = epo_off_c.tolist()
    subj_off_c = op.join(subj + '_off_c')
    subj_130_c = op.join(subj + '_130_c')
    subj_60_c = op.join(subj + '_60_c')
    subj_off_e = op.join(subj + '_off_e')
    subj_130_e = op.join(subj + '_130_e')
    subj_60_e = op.join(subj + '_60_e')
    d = {'subj': [subj_off_c, subj_off_e, subj_130_c, subj_130_e, subj_60_c, subj_60_e],
         'theta': [epo_off_c, epo_off_e, epo_130_c, epo_130_e, epo_60_c, epo_60_e]}
    subj_results = pd.DataFrame(data=d)
    all = all.append(subj_results)
    all = all.explode('theta')

# saving the excel
all.to_excel(file_name)
