
import os
import mne
import os.path as op
from mne import pick_types
from config import eeg_path, subjects

Filter = [1, 48]
reject = dict(eeg=180e-6)
sfreq = 250.  # matlabfiles are already re-sampled!

#define raw brainvision file for INFO, montage file, and "raw", but preprocessed matlab file (DBS filter already applied)
sub_path = op.join(eeg_path, subj)
original_file = os.path.join(sub_path,subj +'_eyelink_'+cond+'.vhdr')
annot_file = os.path.join(sub_path,subj +'_eyelink_'+cond+'.vmrk')
raw_file = os.path.join(sub_path, 'data_rsp_'+ subj +'_'+cond+'_eyelink.mat')
montage_file = os.path.join(eeg_path,'AS-128.bvef')
outfname = op.join(sub_path, subj+'_' + cond+'_raw.fif')

#load original brainvision file to get info
original_data = mne.io.read_raw_brainvision(original_file)
original_info = original_data.info

#read matlab file with brainvision info
raw=mne.io.read_raw_fieldtrip(raw_file, original_info, data_name='data_rsp')
# find events from annotations in original brainvision file - This part is only necessary here if you still need to adjust trial nr in eyelink file
events, event_dict = mne.events_from_annotations(original_data,event_id={'Stimulus/S  2': 2, 'Stimulus/S  3': 3, 'Stimulus/S  9': 9, 'Stimulus/S 99': 99})
# look for pauses in the data collection to annotate these time spans as bad (below when checking for
mne.viz.plot_events(events, event_id=event_dict, sfreq=raw.info['sfreq'], first_samp=raw.first_samp)
# remove accelerometer data - needs to be done in some, but not all data sets
raw.drop_channels(['AccX', 'AccY', 'AccZ'])
# change EKG channel to "official" ECG channel (would be handled as EEG channel otherwise)
raw.set_channel_types({'EKG': 'ecg'})
# set montage
montage = mne.channels.read_custom_montage(montage_file)
raw = raw.set_montage(montage)
# set reference to average
raw.set_eeg_reference('average')

# resample and filter EEG data
picks_eeg = pick_types(raw.info, meg=False, eeg=True, eog=False, emg=False, misc=False, ecg=True, stim=False, exclude='bads')
raw.filter(Filter[0], Filter[1], picks=picks_eeg)

#look for bad channels and add new bad annotation to messy periods
raw.plot_psd(fmax=30, reject_by_annotation=True)
anterior = mne.pick_channels_regexp(raw.ch_names, regexp='A')
raw.plot(order=anterior, n_channels=len(anterior), scalings='auto')
frontal = mne.pick_channels_regexp(raw.ch_names, regexp='F')
fig=raw.plot(order=frontal, n_channels=len(frontal), scalings ='auto')
fig.canvas.key_press_event('a')
central = mne.pick_channels_regexp(raw.ch_names, regexp='C')
fig=raw.plot(order=central, n_channels=len(central), scalings ='auto')
 parietal = mne.pick_channels_regexp(raw.ch_names, regexp='P')
 raw.plot(order=parietal, n_channels=len(parietal), scalings ='auto')
 temporal = mne.pick_channels_regexp(raw.ch_names, regexp='T')
 raw.plot(order=temporal, n_channels=len(temporal), scalings ='auto')
 occipital = mne.pick_channels_regexp(raw.ch_names, regexp='O')
 raw.plot(order=occipital, n_channels=len(occipital), scalings ='auto')

#raw.info['bads'].append('O10')    #add a single channel
#raw.info['bads'].extend(['PO4', 'T7'])  # add a list of channel

# interpolation of bad EEG channels using the spherical spline method
eeg_data = raw.copy().pick_types(meg=False, eeg=True, ecg=True, exclude=[])
eeg_data_interp = eeg_data.copy().interpolate_bads(reset_bads=True)
# check
eeg_data_interp.plot_psd(fmax=30)

raw = eeg_data_interp
# save preprocessed raw file
raw.save(outfname, overwrite=True)

# release memory
del raw
del eeg_data
