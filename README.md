# DBS_f_antisac

TASK

Files required to run the antisaccade task in Matlab with Psychtoolbox are under ./task 
antisac.m is the main scirpt, all others are related to connecting and synchronization of the Eyelink 1000 

ANALYSIS CODE

The R code to break down the raw eye-tracking files (converted into .asc format) into the measures reported in the study 
(latency, error rate, express rate) is under ./eye-tracking
./EEG-analysis contains the Python code required to run all steps of the sensor-level EEG analysis including group-level statistics.
R code for trial-by-trial regressions of band power and antisaccade measures can be found in ./stats


