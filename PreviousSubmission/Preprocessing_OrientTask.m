function Preprocessing_OrientTask(exp)

% Not set-up to process data from the staircasing task

try
    nparts = length(exp.participants);
    nsets = length(exp.setname);
    
    % Replicating event triggers when exp.events is a matrix
    if isempty(exp.epochs) == 1
        exp.epochs = cellstr(num2str(cell2mat(reshape(exp.events,1,size(exp.events,1)*size(exp.events,2)) )'))';
    else
        exp.events = exp.epochs;
    end

    % Is epoch names not specified, use event names
    if isempty(exp.epochs_name) == 1
        exp.epochs_name = exp.event_names;
    else
        exp.event_names = exp.epochs_name;
    end
    
    
    for i_set = 1:nsets
        
        sprintf(exp.setname{i_set})
        
        % if folder doesn't exist yet, create one
        if ~exist([exp.dataLocation  '\ProcessData\' exp.setname{i_set} '\Segments\'])
            mkdir([exp.dataLocation  '\ProcessData\' exp.setname{i_set} '\Segments\']);
        end
        
        %initialize EEGLAB
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %subject numbers to analyze
        
        nevents = length(exp.events(i_set,:));
        
        %% Load data and channel locations
        for i_part = 1:nparts
%         for i_part = 12    
            
            part_name = ['subj' exp.participants{i_part}]; %this is cuz subject ids are in the form 'subj1' rather than '001' in orientation wheel
            
            if strcmp('on',exp.epoch) == 1
                
                sprintf(['Participant ' num2str(exp.participants{i_part})])
                
                %% Load a data file
                if strcmpi(part_name,'subj13') %subject's EEG file has different name
                    EEG = pop_loadbv(exp.pathname, [part_name '_staircase.vhdr']); %orientation wheel data only
                else
                    EEG = pop_loadbv(exp.pathname, [part_name '_orient.vhdr']); %orientation wheel data only
                end
                
                %% Load channel information
                EEG = pop_chanedit(EEG, 'load',{exp.electrode_locs 'filetype' 'autodetect'});
               
                %% Arithmetically re-reference to linked mastoid (M1 + M2)/2
                % only for brain electrodes (exclude mastoids & EOGs)
                for ii = exp.brainelecs(1):length(exp.brainelecs)
                    EEG.data(ii,:) = (EEG.data(ii,:)-((EEG.data(exp.refelec,:))*.5));
                end
                clear ii
                
                %% Filter the data
                if strcmpi(exp.filter,'on')
%                    EEG = pop_eegfilt( EEG, 0, 30, [], 0); %with low pass of 30
                   EEG = pop_eegfiltnew(EEG, exp.locutoff, exp.hicutoff); % filter function
                end
                
                %% Change markers so they can be used by the gratton_emcp script
                allevents = length(EEG.event);
                for i_event = 2:allevents %skip the first
                    EEG.event(i_event).type = num2str(str2num(EEG.event(i_event).type(2:end)));
                end

                %% The triggers are early
                [EEG] = VpixxEarlyTriggerFix(EEG);
                
                %% Extract epochs of data time locked to event
                %Extract data time locked to targets and remove all other events
                EEG = pop_epoch(EEG, exp.epochs, exp.epochslims, 'newname', [part_name '_epochs'], 'epochinfo', 'yes');
                %subtract baseline
                EEG = pop_rmbase(EEG, exp.epochbaseline);
  
                %% Get behavior data and add to EEG structure
                [error_deg] = getBEHdata_OrientTask(part_name,exp);
                % add error deg to epoch structure
                if length(EEG.epoch) == length(error_deg)... %make sure right BEH file
                        || strcmpi(part_name,'subj6') || strcmpi(part_name,'subj8') %or subject 6 and 8
                   EEG.error_deg = error_deg;
                end
                
                %% Reject practice trials from data
                % only specific to subject 8 % only when using these specific settings
                if (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byTargets'))...
                    || (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byCatchTargets'))
                    rej_practice = zeros(1,length(EEG.epoch));
                    rej_practice(1,1:18) = 1; %mark the first 18 trials for removal
                    EEG = pop_rejepoch(EEG, rej_practice, 0);
                end
                clear rej_practice
                

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %% Artifact Rejection, EMCP Correction, then 2nd Rejection
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 1, trials with range >exp.preocularthresh uV
                if isempty(exp.preocularthresh) == 0
                    rejtrial = struct([]);
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],exp.preocularthresh(1),exp.preocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    % `````````````````````````````````````````````````````   
                    rejtrial(i_set,1).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % EMCP occular correction
                temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
                EEG = gratton_emcp(EEG, exp.selection_cards, {'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
                EEG.emcp.table %this prints out the regression coefficients
                EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Baseline again since EMCP changed it
                EEG = pop_rmbase(EEG,exp.epochbaseline);
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 2, trials with range >exp.postocularthresh uV
                if isempty(exp.postocularthresh) == 0
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],exp.postocularthresh(1),exp.postocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    rejtrial(i_set,2).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
                %% Additional rejection of trials (reviewed visually) 
                % Need to set breakpoint at below function to be able to 
                % make changes to this code without exiting this function
                    % Use pop_eegplot(EEG) to visualize activity
                EEG = reject_trials_inspect(EEG,part_name,exp,rejtrial,i_set);
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                %% Replace the stored data with this new set
                tempEEG = EEG;               
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                 
                %% Select individual events
                for i_event = 1:nevents
                    EEG = pop_selectevent( tempEEG, 'type', exp.events{i_set,i_event}, 'deleteevents','on','deleteepochs','on','invertepochs','off');
                    EEG = pop_editset(EEG, 'setname', [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}] );
                    EEG = pop_saveset( EEG, 'filename',[part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.set'],'filepath',[exp.pathname '\' exp.setname{i_set} '\Segments\']);
                end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                 
            end %create epoch loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
            %% Time-Frequency Data
            if strcmp('on',exp.tf) == 1 || strcmp('on',exp.singletrials) == 1

                if ~exist([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\'])
                    mkdir([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\']);
                end

                
                for i_event = 1:nevents
                    
                    if strcmp('on',exp.epoch) == 0 %loading previous epochs if not created this session
                        filename = [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
                        EEG = pop_loadset('filename',[filename '.set'],'filepath',...
                            [exp.dataLocation '\ProcessData\' exp.setname{i_set} '\Segments\']);
                    end
                    
                    if isempty(exp.tfelecs) %if TF electrodes not specified, same as exp.brainelecs
                       exp.tfelecs = exp.brainelecs; 
                    end
                    
                    for i_tf = 1:length(exp.tfelecs)
                        i_chan = exp.tfelecs(i_tf);
                        EEG = eeg_checkset(EEG);
                        [ersp(i_chan,:,:),itc(i_chan,:,:),powbase,times,freqs,dum1,dum2,all_ersp(i_chan).trials] =...
                            pop_newtimef(EEG, 1, i_chan, exp.epochslims*1000, exp.cycles, ...
                            'topovec', i_chan, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo,...
                            'baseline', exp.erspbaseline, 'freqs', exp.freqrange, 'freqscale', 'linear', ...
                            'padratio', exp.padratio, 'plotphase','off','plotitc','off','plotersp','off',...
                            'winsize',exp.winsize,'timesout',exp.timesout);
                    end
                    clear i_chan i_tf

                    if strcmp('on',exp.tf) == 1 %if TF was done already, do not save
                        if exp.cycles(1) > 0 %for wavelet
                            save([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\Wav_' part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'ersp','itc','times','freqs','powbase','exp')
                        else %for FFT
                            save([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\TimeFrequency\TF_' part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'ersp','itc','times','freqs','powbase','exp')
                        end
                    end
                        
                    
                     % Save single trial data
                    if strcmp('on',exp.singletrials) == 1
                        
                        % Create folder for single trial data
                        if ~exist([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\'],'dir')
                            mkdir([exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\']);
                        end
                        
                        % File path name
                        Filepath_Trials = [exp.dataLocation '\ProcessData\' exp.setname{i_set} '\SingleTrials\' part_name '\'];
                        
                        % Save single trial data from the selected electrodes
                        for zzz = 1:length(exp.singletrialselecs)
                            i_chan = exp.singletrialselecs(zzz);
                            elec_all_ersp = all_ersp(i_chan).trials;
                            if exp.cycles(1) > 0 %for wavelet
                                save([Filepath_Trials exp.singtrlelec_name{zzz} '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '_Wav.mat'],...
                                'elec_all_ersp','times','freqs','powbase','exp')
                            else %for FFT
                                save([Filepath_Trials exp.singtrlelec_name{zzz} '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],...
                                'elec_all_ersp','times','freqs','powbase','exp')
                            end
                        end
                        clear i_chan elec_all_ersp zzz
                    end
                    clear Filepath_Trials

                end
                eeglab redraw
                
            end
            clear part_name
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             
        end %i_part loop end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::         
    end %i_set loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    
% !!!!!!!!!!!!!!!!    
catch ME % !!!!!!!
    save('dump') 
    throw(ME)
end % !!!!!!!!!!!!
% !!!!!!!!!!!!!!!! 


% ####################
end % end function ###
% ####################
