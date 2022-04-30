# PD_DBS_Antisac

High- and Low-Frequency Deep Brain Stimulation in the Subthalamic Nucleus differentially modulate Response Inhibition and Action Selection in Parkinson’s Disease

Josefine Waldthaler, Alexander Sperlich, Aylin König, Charlotte Stüssel, Frank Bremmer, Lars Timmermann & David Pedrosa 

This repo contains the code required to reproduce the task and most of the analysis presented in this study: 


TASK

Files required to run the antisaccade task in Matlab with Psychtoolbox are under ./task 
./task/antisac.m is the main scirpt, the others are related to connecting and synchronizing of the Eyelink 1000 

ANALYSIS CODE

The R code to break down the raw eye-tracking files (converted into .asc format) into the measures reported in the study 
(latency, error rate, express rate) is under ./eye-tracking
./EEG-analysis contains the Matlab code for the first preprocessing steps (Hampel filtering to remove DBS artefacts) and the Python code required to run all steps of the sensor-level EEG analysis including group-level statistics.
R code for trial-by-trial regressions of band power and antisaccade measures can be found in ./stats

DATA

The raw datasets cannot be made publicly available due to the German and European data privacy laws, but can be requested from the corresponding author in an anonymized form for replication purposes. 
