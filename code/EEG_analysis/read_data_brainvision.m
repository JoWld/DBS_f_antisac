function read_data_brainvision(wdir, subj, experiment)

%   This function reads (raw) data and performs a resampling with the
%   fieldtrip toolbox; besides, artefact correction for DBS signals
%   is performed on the raw signal.

addpath /Users/jowld/Documents/EEG_DBS_Inhib/matlab/fieldtrip
ft_defaults

ROOTDIR = '/Users/jowld/Documents/EEG_DBS_Inhib/EEG_data/matlab';
wdir = '/Users/jowld/Documents/EEG_DBS_Inhib/EEG_data';     
cd(fullfile(wdir))

subj            = "52_MKA_1308"

%% General settings
rspl_freq       = 250;                                                      % frequency at which the data will be resampled
conds           = {'off', '130'}; %'60', '180'
nchans          = 132;
debug           = 0;
experiment      = "eyelink";

outdir = fullfile(wdir, 'subj');                                        % directory in which data will be saved

% Start extracting and resampling data
for s = 1:numel(subj)
    %fprintf('\n\tprocessing subj: %s\n', num2str(subj{s}))
    dir_rawdata = fullfile(wdir, subj);                                  % folder in which raw data is stored
    % file_prefix = subjdetails(ROOTDIR, subj);
    
    for c = 1:numel(conds)
        fprintf('\t\t ... condition: %s', num2str(upper(conds{c})))
        filename_save = strcat('data_rsp_',subj, ...
            conds{c}, sprintf('_%s', experiment), '.mat');                  % filename under which data will be saved in the (outdir) director
         
        filename2load = fullfile(dir_rawdata, strcat(subj, ...
                sprintf('_%s_', experiment), conds{c}, '.eeg'));            % filename, so that data may be loaded
            
%             if ~exist(filename2load, 'file')
%                 fprintf('\t\tproblem with reading data from subj: %s, cond: %s. Please select file manually \n', s, conds{c});
%                 
%                 cd(dir_rawdata)                                             % change directory and select file manually
%                 [file,path] = uigetfile('*.eeg');
%                 if isequal(file,0)
%                     fprintf('No file selected, skipping to next subject ...\n');
%                     continue
% %                 else
%                     filename2load = file;
%                     fprintf('\n\t\tThe selected file for the %s condition is: %s\n', conds{c}, fullfile(path,file));
%                 end
%                 cd(wdir
       
            %   Reads data from original brainvision-files taking
            %   specific channels into account
            cfg = [];                                                       % cfg is used in the for-loop for reading and downsampling data channelwise
            cfg.resamplefs = rspl_freq;
            cfg.detrend    = 'no';
            cfg.feedback   = 'no';
            
            singlechan = cell(1, nchans);  %cell(1,4)%                      % pre-allocate space
            for i=1:nchans % for loop is necessary because of size if the data over 2GB
                cfg_temp = [];                  
                cfg_temp.dataset    = filename2load;       % can bei either string (channel name) or number
                cfg_temp.channel    = i;                    % reads data from filename as defined before
                temp                = ft_preprocessing(cfg_temp);
%                 fprintf('\n\t\t\t... processing %s {%s of %s channels = %.1f%% done} \n', ...
%                     temp.label{1}, num2str(i-1), num2str(nchans), (i-1)/nchans*100)

               % remove DBS artfact
                if ~strcmp(conds{c}, 'off')
                    temp_filtered = ...
                        DBSartefacts_removal(temp, 'simple_filter', ...
                        str2double(conds{c}));
                    temp_processed = ...
                        DBSartefacts_removal(temp_filtered, ...
                        'hampel_identifier', str2double(conds{c}));
                    
                    if debug == 1
                        plot_removeDBSartefact(temp, temp_processed, ...
                            str2double(conds{c}), 100)
                    end
                else
                    temp_processed = temp;
                end
                singlechan{i}       = ft_resampledata(cfg, temp_processed);
            end
            fprintf('\n\t\tprocessed all %s channels\n', num2str(nchans))
            
            % Append data to one file/
            cfg = [];
            data_rsp = ft_appenddata(cfg, singlechan{:});
            
            if debug == 1
                cfg = [];
                cfg.viewmode = 'vertical';
                ft_databrowser(cfg, data_rsp)
            end
            
            % Saves data to pre-specified folder
            save(fullfile(outdir, filename_save), 'data_rsp', '-v7.3');
            

        end
    end
end