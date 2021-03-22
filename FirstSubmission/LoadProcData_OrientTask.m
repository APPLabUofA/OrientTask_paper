% This code is for loading processed EEG data. Settings below are to
% specify what type of data
% Use processing settings to specify which processed data set to load. 

% -------------------------------------------------------------------------
%% Load processing settings
% load('filt_byCatchTargets_v2_wav_Settings.mat'); %this is to load target-aligned catch trials
load('filt_byTargets_v4_Settings.mat'); %this loads all processed target-aligned trials except catch trials
% -------------------------------------------------------------------------
%% Do not use if loading catch trials
anal.rm_rejtrial_BEH = 'off'; % return response error w/rejected trials removed & model fit 
% -------------------------------------------------------------------------
%% Select data by type of TF analysis: 'FFT' or 'wavelet' 
anal.method = 'wavelet'; %published processed data is only this type
% anal.method = 'FFT'; 
% -------------------------------------------------------------------------
%% Specify type of data to load
% (epoched segements, single-trial,time-frequency, etc...) 
anal.tf = 'off'; % if loading TF data
anal.singletrials = 'on'; % if loading single trial data
anal.segments = 'off'; % if loading epochs
anal.tfelecs = exp.brainelecs; %which TF electrodes if only loading a few
anal.singletrialselecs = exp.singletrialselecs; %which single trial electrodes if loading a few
% -------------------------------------------------------------------------


% #########################################################################
% #########################################################################

nparts = length(exp.participants); %number of subjects
nsets = length(exp.setname); %will always be 1 in this study

% -------------------------------------------------------------------------
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% -------------------------------------------------------------------------



% #########################################################################
% #########################################################################
%% Load the data
% The main loop loops through sets, then participants, then events.
for i_set = 1:nsets
    exp.setname(i_set)
    tic
    
%     for i_part = nparts
    for i_part = 1:nparts
        sprintf(['Loading Participant ' num2str(exp.participants{i_part}) '...' ])
        
        % number of events of interest (not really needed here except when
        % there are muliple sets)
        nevents = length(exp.events(i_set,:));
        
        part_name = ['subj' exp.participants{i_part}]; %this is cuz subject ids are in the form 'subj1' rather than '001'
        
        for i_event = 1:nevents
            
            filename = [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
            
            % --------------------------------------------------------------------------------------------------------------------
            % Load the Time frequency data, if needed.
            if strcmpi('on',anal.tf) == 1 % only load these variables if we are loading time-frequency data
                
                if strcmpi('FFT',anal.method) == 1 % load FFT data
                    %The variable ersp will be a 6D variable: (participants x sets x events x electrodes x frequencies x timepoints).
                    ersp(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'ersp'));
                    itc(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'itc'));

                    if i_part == 1 && i_set == 1 && i_event == 1 %load time and freq data
                        times = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'times'));
                        freqs = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'freqs'));
                    end
                    
                elseif strcmpi('wavelet',anal.method) == 1 % load wavelet data
                    %The variable ersp will be a 6D variable: (participants x sets x events x electrodes x frequencies x timepoints).
                    ersp(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'ersp'));
                    itc(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'itc'));

                    if i_part == 1 && i_set == 1 && i_event == 1 %load time and freq data
                        times = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'times'));
                        freqs = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'freqs'));
                    end
                end
                
            end
            % --------------------------------------------------------------------------------------------------------------------
    
            % --------------------------------------------------------------------------------------------------------------------
            % Load the EEGLAB datasets, if needed.
            if strcmpi('on',anal.segments) == 1 || strcmp('on',anal.singletrials) == 1 % only load these variables if we are loading either ERP or single trial data
                try
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.dataLocation '\ProcessData\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                catch
                    WaitSecs(.5)
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.dataLocation '\ProcessData\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            elseif strcmpi('on',anal.tf) == 1 % if we are loading time-frequency data only, then we just need one of these.
                if i_part == 1 && i_set == 1 && i_event == 1
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.dataLocation '\ProcessData\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            end
            % --------------------------------------------------------------------------------------------------------------------
            
            % --------------------------------------------------------------------------------------------------------------------
            % Load the Single Trial complex values, if needed
            if strcmpi('on',anal.singletrials) == 1 % only load these variables if we are loading single trial data
                
%                 ntrigs = length(exp.events{i_set});
                
                % Loads the time values and freqs (freqs will be different depending on TF method)
                if strcmpi('FFT',anal.method) == 1 % load FFT data
                    if i_part == 1 && i_set == 1 && i_event == 1
                        times = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'times'));
                        freqs = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'TF_' filename '.mat'],'freqs'));
                    end
                elseif strcmpi('wavelet',anal.method) == 1 % load wavelet data
                    if i_part == 1 && i_set == 1 && i_event == 1
                        times = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'times'));
                        freqs = struct2array(load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\' 'Wav_' filename '.mat'],'freqs'));
                    end
                end
                

                % Load single trial data for each electrode
                for ii = 1:length(exp.singletrialselecs)
                    i_chan = exp.singletrialselecs(ii);
                    
                    % all_ersp is (participant x electrode).trials(freq x time x trial)
                    if strcmpi('FFT',anal.method) == 1 % load FFT data
                        try %Unfortunately, this load procedure can break sometimes in a non-reproducible way. So if an error happens here, we wait half a second and try again.
                            channeldata = load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
                            all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        catch
                            WaitSecs(.5)
                            channeldata = load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
                            all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        end
                    elseif strcmpi('wavelet',anal.method) == 1 % load wavelet data
                         try %Unfortunately, this load procedure can break sometimes in a non-reproducible way. So if an error happens here, we wait half a second and try again.
                            channeldata = load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '_Wav.mat'],'elec_all_ersp');
                            all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        catch
                            WaitSecs(.5)
                            channeldata = load([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '_Wav.mat'],'elec_all_ersp');
                            all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        end
                    end
                    clear channeldata i_chan
                    
                 end
            % --------------------------------------------------------------------------------------------------------------------
            end
        end
    end
    toc
end
clear i_event i_set i_part filename nevents nparts nsets ii part_name


eeglab redraw


% #########################################################################
% #########################################################################
%% Return BEH data corrected for rejected trials
% Get degree response errors on trials with the trials rejected during the
% EEG pre-processing removed & model fit parameters
if strcmpi('on',anal.rm_rejtrial_BEH) == 1 &&...
    strcmpi('filt_byCatchTargets_v2_wav',exp.setname) == 0 %do not use on catch trials
    [resp_errdeg, model_out] = getRespErr_ModelFits(exp,ALLEEG);
    
    %save processed behavior data and epoch information
    save([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'],'resp_errdeg','model_out','freqs','times')
end


% #########################################################################
% -------------------------------------------------------------------------
% #########################################################################

clear anal

% #########################################################################
% -------------------------------------------------------------------------
% #########################################################################


