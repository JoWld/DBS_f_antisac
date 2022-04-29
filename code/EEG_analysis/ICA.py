#   This script is part (3) of the EEG analysis of an antisaccade task and 
#   runs an ICA.

#  Copyright (C) April 2020, last modified April 2022
#   J.Waldthaler, A. Sperlich, D. Pedrosa
#   University Hospital of Gießen and Marburg
#
#   This software may be used, copied, or redistributed as long as it is
#   not sold and this copyright notice is reproduced on each copy made.
#   This routine is provided as is without any express or implied
#   warranties whatsoever.


from mne.preprocessing import (ICA, create_ecg_epochs)
import numpy as np
from os import mkdir
from config import subjects, eeg_path, sub_path

# loop if wanted
#for subj in subjects_and_dates:
    print('Now processing subject ' + subj)
    # define paths
    sub_path = op.join(eeg_path, subj)
    ica_path = op.join(eeg_path, subj, 'ica')
    raw_file = os.path.join(sub_path, subj+'_' + cond+'_raw.fif')
    out_icaFname = op.join(ica_path, cond+'_comp-ica.fif')  # Name of ICA components (saved for bookkeeping)
    outfname = op.join(sub_path, subj +'_' + cond+ '_ica-raw.fif')

    if op.exists(outfname):
        print('Output ' + outfname + ' already exists')
    # make dirs if they do not exist
    if not op.exists(sub_path):
        mkdir(sub_path)
    if not op.exists(ica_path):
        mkdir(ica_path)

    # load data
    raw = mne.io.read_raw(raw_file, preload=True)
    picks_eeg = pick_types(raw.info, meg=False, eeg=True, eog=False, emg=False, misc=False, ecg=False,
                       stim=False, exclude='bads')
    # run ICA
    ica = ICA(n_components=0.95, method='fastica', random_state=0)
    ica.fit(raw, picks=picks_eeg, decim=3, verbose=True, reject_by_annotation=True)

    # plot and save ICA components
    ica_fig = ica.plot_components()
    [fig.savefig(op.join(ica_path, 'ICA_allComp' + str(i) + '.png')) for i, fig in enumerate(ica_fig)]
    ica.save(out_icaFname)
    print('ICA comp saved as ' + out_icaFname)
    # inspect ICA components and find blinks / eye movements
    ica.plot_sources(raw)

    # choose components to exclude
    #ica.exclude = [0,3,5]

    # remove ECG component
    n_max_ecg = 3
    picks_ecg = mne.pick_types(raw.info, meg=False, eeg=False, eog=False, emg=False, misc=False, ecg=True, stim=False, exclude='bads')
    raw.notch_filter(50,  picks=picks_ecg)  # Remove residual 50Hz line noise, sometimes causes trouble

    # Find ECG artifacts
    ecg_epochs = create_ecg_epochs(raw, ch_name='EKG', tmin=-.5, tmax=.5)  # , picks=picks)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, method='ctps')
    # Update reject info
    ica.exclude += ecg_inds[:n_max_ecg]
    # Plot ECG ICs for inspection
    ecg_scores_fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, title='Component score (ecg)', show=True)
    ecg_scores_fig.savefig(op.join(ica_path, 'ICA_ecg_comp_score.png'))
    #plt.close()

    if ecg_inds:
        show_picks = np.abs(ecg_scores).argsort()[::-1][:5]
        ecg_comp_fig = ica.plot_components(ecg_inds, title='ecg comp', colorbar=True, show=False)
        ecg_comp_fig.savefig(op.join(ica_path, 'ICA_ecg_comp_topo.png'))
        #plt.close()
        # estimate average artifact
        ecg_evoked = ecg_epochs.average()
        # plot ECG sources + selection
        ecg_evo_fig1 = ica.plot_overlay(ecg_evoked, exclude=ecg_inds, show=False)
        ecg_evo_fig1.savefig(op.join(ica_path, 'ICA_ecg_overlay.png'))

    # apply the ICA solution to raw and save file
    raw_ica = ica.apply(raw)
    raw_ica.save(outfname, overwrite=True)
    print('----------- FINISHED ' + subj + ' -----------------')

# release some memory
del raw_ica
del ica
