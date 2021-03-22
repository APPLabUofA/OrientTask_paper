% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%%                       Make ERPs by SD Parameter
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% **Need to have loaded the data using LoadProcData_OrientTask.m**
% Code below is for regular trials (not catch trials) that has been aligned
% to target onset

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('filt_byTargets_v4_Settings.mat'); %setting for target-aligned trials except catch trials

% /////////////////////////////////////////////////////////////////////////
%% Location to save ERP results
saveLocation = [exp.dataLocation '\ProcessData\']; % set save directory of data set

% /////////////////////////////////////////////////////////////////////////
%% Set electrodes to analyze if only doing a few 
% (default is all electrodes)
elect_erp = 1:EEG.nbchan;

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////
%% Creates ERPs based on SD parameter from model
errlims = cell(1,length(exp.participants));    %pre-allocate
errs_x = cell(1,length(exp.participants));    %pre-allocate
errs_n = cell(1,length(exp.participants));    %pre-allocate
erp_out_x = NaN(length(exp.participants),length(elect_erp),length(EEG.times)); %pre-allocate
erp_out_n = NaN(length(exp.participants),length(elect_erp),length(EEG.times)); %pre-allocate

% create separate ERPs for large and small errors
for i_part = 1:length(exp.participants)
    
    % Get upper and lower limits based on model fit
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    % Calculate ERP
    for ii = 1:length(elect_erp)
        i_elect = elect_erp(ii); %for doing only a selection of electrodes
        
        % Get trials with small errors
        erp_out_x(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
            [find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] ),3));
        %save small errors
        errs_x{i_part} = resp_errdeg{i_part}(resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75));
        
        % Get trials with large errors
        erp_out_n(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,...
            [find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] ),3));
        %save large errors
        errs_n{i_part} = resp_errdeg{i_part}([find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))]);
    end
    clear ii i_elect
end
clear i_part

% /////////////////////////////////////////////////////////////////////////
% Save data
erp_time = EEG.times; %save time variable
chan_locs = EEG.chanlocs; %save channel locations
save([saveLocation 'erp_out_byTarget.mat'],'errs_x','errs_n','erp_out_x','erp_out_n',...
    'model_out','resp_errdeg','erp_time','chan_locs')

% /////////////////////////////////////////////////////////////////////////
% Clear workspace
ccc
% /////////////////////////////////////////////////////////////////////////





% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%%                    Make ERPs for Catch Trials
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% **Need to have loaded the data using LoadProcData_OrientTask.m**
% Code below is for catch trials that has been aligned to target onset

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('filt_byCatchTargets_v2_wav_Settings.mat'); %setting for target-aligned trials except catch trials

% /////////////////////////////////////////////////////////////////////////
%% Location to save ERP results
saveLocation = [exp.dataLocation '\ProcessData\']; % set save directory of data set

% /////////////////////////////////////////////////////////////////////////
%% Set electrodes to analyze if only doing a few 
% (default is all electrodes)
elect_erp = 1:EEG.nbchan;

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////
%% Creates ERPs averaged over catch trials

erp_out = NaN(length(exp.participants),length(elect_erp),length(EEG.times)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(elect_erp)
        i_elect = elect_erp(ii); %for doing only a selection of electrodes
        erp_out(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,:),3));
        clear i_elect
    end
    clear ii
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
% Save data
erp_time = EEG.times; %save time variable
chan_locs = EEG.chanlocs; %save channel locations
save([saveLocation 'erp_out_catch_trials.mat'],'erp_out','erp_time','chan_locs')

% /////////////////////////////////////////////////////////////////////////
% Clear workspace
ccc
% /////////////////////////////////////////////////////////////////////////





