%% ANALYSIS WRAPPER

% ==== TRIGGERS ====
%
% -- Targets Present --
% Trial start (fixation): lag (1,2,3,4)
% Entrainers: 61-68
% Target: 20+lag (21,22,23,24)
% Mask: 90+lag (91,92,93,94)
% Response screen: 40+lag (41,42,43,44)  
% Response: 80+lag (81,82,83,84)
%
% -- Targets Not Present --
% Trial start (fixation): 10+lag (11,12,13,14)
% Entrainers: 61-68
% Target: 30+lag (31,32,33,34)
% Mask: 95+lag (96,97,98,99)
% Response screen: 50+lag (51,52,53,54)%
% Response: 70+lag (71,72,73,74)



%clear and close everything
ccc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        >>>>> Description of the saved dataset settings <<<<<
% 
% 
% -> filt_byTargets_v4: winsize is 256, no ERSP baseline, epoched to targets; 
%    Epoch limit [-1.5 1.5]. ERP baseline [-200 0]. Filter on [0.1 50]. 
%    Cycles [2 0.6] (2 cycles at lowest freq & 12 at highest). Freq range [2 40]. 
%    Freq increase in steps of 1.027 Hz. Timesout = 300. Padratio = 4. 
%    Subjects removed (subj 9 (7) & subj 22 (20))
% 
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> filt_byCatchTargets_v2_wav: winsize is 256, no ERSP baseline. Epoched to 
%    targets on catch trials; Epoch limit [-1.5 1.5]. ERP baseline [-200 0]. 
%    Filter on [0.1 50]. Cycles [2 0.8] (2 cycles at lowest freq & 8 at 
%    highest). Freq range [2 40]. Timesout = 300. Padratio = 4.  
%    Subjects removed (subj 9 (7) & subj 22 (20))
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set file locations
exp.dataLocation = [pwd '\Data_Files']; % set directory of data set (default: pwd)

%% Load data settings
exp.name = 'OrientWheel';
exp.conds = ''; %conds is usually used for comparing the same type of trials
                %under different conditions (e.g., stimulation vs sham)
% exp.pathname = 'M:\Data\OrientationWheel\EEG\'; %path of EEG data
exp.pathname = [exp.dataLocation '\RawData\EEG\']; %path of EEG data
% exp.setname = {'filt_byTargets_v4'}; % name each epoched set
exp.setname = {'filt_byCatchTargets_v2_wav'}; % name each epoched set
% note: the meaning of set here is to identify the type of analysis done.
%       set is usually used to identify different trial types (e.g., standards
%       vs targets) within the same experimental condition.

% List of participants' ids
% note: removed subj 23 dues to excessive movement
% exp.participants = {'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'24';...
%     '25';'26';'27';'29';'30'};
exp.participants = {'3';'4';'5';'6';'7';'8';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'24';...
    '25';'26';'27';'29';'30'}; %subj 9 & subj 20 removed

%% Blink Correction
% the Blink Correction wants dissimilar events (different erps) seperated by 
% commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
% e.g., exp.selection_cards = {'11 21','13 23'}; %must be list == length(exp.setname)
%%%indicates where you want to center your data (where time zero is)
%must be list == length(exp.setname)
% exp.selection_cards = {'81 82 83 84'}; %response aligned 
% exp.selection_cards = {'21 22 23 24'}; %target
exp.selection_cards = {'31 32 33 34'}; %catch trials
% exp.selection_cards = {'41 42 43 44'}; %response screen aligned 

%% Artifact rejection. 
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold 
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
% exp.postocularthresh = [ ]; %Second happens after. Leave blank [] to skip
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip

%% Events and event labels
%events are what happen within each trial (e.g., fixation, target, response, etc...) 
%%%for each condition (lag 1-4 for example), numbers correspond to
%%%triggers that will be kept for each condition. All other triggers will
%%%be removed
%can be list or matrix (sets x events) 
% exp.events = {[81,82,83,84]};%response aligned 
% exp.events = {[21,22,23,24]};%target-aligned 
exp.events = {[31,32,33,34]};%target-aligned catch trials 
% exp.events = {[41,42,43,44]};%response screen aligned

% exp.event_names = {'Wheel','Wheel','Wheel','Wheel'}; %must be list or matrix (sets x events)
% exp.event_names = {'Target','Target','Target','Target'}; %must be list or matrix (sets x events)
exp.event_names = {'cTarget','cTarget','cTarget','cTarget'}; %must be list or matrix (sets x events)
% exp.event_names = {'Resp','Resp','Resp','Resp'}; %must be list or matrix (sets x events)

% exp.suffix = {'byTarg'};
exp.suffix = {'bycTarg'}; %if target was present 

%% Electrode location
%Where are your electrodes? (.ced file)
exp.electrode_locs = [pwd '\EOG-electrode-locs-32_orientwheel.ced'];
% electrode information
exp.electrode = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
exp.elec_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

%% Re-referencing the data
exp.refelec = 1; %which electrode do you want to re-reference to?
exp.brainelecs = [2:32]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Filter the data?
exp.filter = 'on'; %filter all files
exp.hicutoff = 50; %higher edge of the frequency pass band (Hz)
exp.locutoff = 0.1; %lower edge of the frequency pass band (Hz)

%% FFT/Wavelet settings
% How long is your window going to be? (Longer window == BETTER frequency 
% resolution & WORSE time resolution)
exp.winsize = 512; %use numbers that are 2^x, e.g. 2^10 == 1024ms

% Baseline will be subtracted from the power variable. It is relative to 
% your window size. Can use just NaN for no baseline
%e.g., [-200 0] will use [-200-exp.winsize/2 0-exp.winsize/2]; 
exp.erspbaseline = NaN;
% exp.erspbaseline = [-400 -200];

% Instead of choosing a windowsize, you can choose a number of cycles per 
% frequency for standard wavelet analysis: usually [3] for better temporal
% precision or [6] for better frequency precision.
% If [wavecycles factor], wavelet cycles increase with frequency beginning 
% at wavecyles. See "help popnewtimef"
% exp.cycles = [0]; %leave it at 0 to use FFT
exp.cycles = [2 0.7]; %number of cycles 

% Choose number of output times
exp.timesout = 300; %200 is usually used

% Set sampling factor for frequencies. 
% when exp.cycles==0, frequency spacing is (low_freq/padratio). For wavelet,
% multiplies the # of output freqs by dividing their spacing (2 is default).
% higher values give smooth looking plot but at computational cost (16 is
% very high)
exp.padratio = 4;

% What frequencies to consider?
% exp.freqrange = [1 40]; 
exp.freqrange = [exp.cycles(1) 40]; %when doing wavelet

%% Epoching the data
exp.epoch = 'on'; %on to epoch data; off to load previous data
%%%indicates where you want to center your data (where time zero is)
exp.epochs = {}; %must be list == length(exp.setname)
exp.epochs_name = {};
exp.epochslims = [-1.5 1.5]; %in seconds; epoched trigger is 0 e.g. [-1 2]
exp.epochbaseline = [-200 0]; %remove the baseline for each epoched set, in ms. e.g. [-200 0] 


%% Time-Frequency settings
%Do you want to run time-frequency analyses? (on/off)
exp.tf = 'on';
%Do you want to use all the electrodes or just a few? Leave blank [] for 
% all (will use same as exp.brainelecs)
exp.tfelecs = [];

%Do you want to save the single-trial data? (on/off) (Memory intensive!!!)
exp.singletrials = 'on';
%Saving the single trial data is memory intensive. Just use the electrodes
% you need. 
exp.singletrialselecs = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
exp.singtrlelec_name = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


%//////////////////////////////////////////////////////////////////////////
%% Save your pipeline settings
% The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
exp.settings = char(exp.setname); %name settings
% `````````````````````````````````````````````````````````````````````````
% Saving will help you remember what settings were used in each dataset
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. 
%//////////////////////////////////////////////////////////////////////////


%% Run preprocessing code
tic %start timer
Preprocessing_OrientTask(exp)
toc %end timer

