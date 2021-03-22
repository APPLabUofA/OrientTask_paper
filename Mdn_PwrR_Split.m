% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% The code below can be run on other processes data (like the normalized
% EEG data). Only the raw log power data was used in the manuscript.
% 
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 
currentfolder = pwd; %to return to current folder after loaded data

% /////////////////////////////////////////////////////////////////////////
%% Load Processed EEG Data
all_erspR = struct2cell(load([saveLocation 'all_ersp_R_v4.mat'],'all_erspR'));  %gets loaded as a struct
all_erspR = all_erspR{1};
chanlocs = struct2cell(load([saveLocation 'all_ersp_R_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;



% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%%             Separate trials based on median band power
% /////////////////////////////////////////////////////////////////////////
% #########################################################################
% This takes a very very long time so be prepared to leave your computer
% running for a few days
errdeg_Hpwr_trl = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
errdeg_Lpwr_trl = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
g_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
g_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
sd_out_Hpwr = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
sd_out_Lpwr = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        tmp_pwr = all_erspR{i_part,i_elect}.trials;
        for i_t = 1:length(times)
            for i_f = 1:length(freqs)

                % Get trials by median power split
                mdn_pwr = nanmedian(tmp_pwr(i_f,i_t,:)); %get single subject's ersp
                
                errdeg_Hpwr_trl{i_part,i_elect}(i_f,i_t,:) = resp_errdeg{i_part}(tmp_pwr(i_f,i_t,:)>=mdn_pwr);
                errdeg_Lpwr_trl{i_part,i_elect}(i_f,i_t,:) = resp_errdeg{i_part}(tmp_pwr(i_f,i_t,:)<mdn_pwr);
               
                if i_part == 21 || i_part == 22 %subjects had bias in mean errors
                    model_out_tmp = MLE(errdeg_Hpwr_trl{i_part,i_elect}(i_f,i_t,:),model2); %fits without plotting
                    model_out_Hpwr_trl = model_out_tmp(2:3);
                    clear model_out_tmp
                    model_out_tmp = MLE(errdeg_Lpwr_trl{i_part,i_elect}(i_f,i_t,:),model2); %fits without plotting
                    model_out_Lpwr_trl = model_out_tmp(2:3);
                    clear model_out_tmp
                else
                    model_out_Hpwr_trl = MLE(errdeg_Hpwr_trl{i_part,i_elect}(i_f,i_t,:),model); %fits without plotting
                    model_out_Lpwr_trl = MLE(errdeg_Lpwr_trl{i_part,i_elect}(i_f,i_t,:),model); %fits without plotting
                end
        
                g_out_Hpwr(i_part,i_elect,i_f,i_t) = model_out_Hpwr_trl(1);
                g_out_Lpwr(i_part,i_elect,i_f,i_t) = model_out_Lpwr_trl(1);
                sd_out_Hpwr(i_part,i_elect,i_f,i_t) = model_out_Hpwr_trl(2);
                sd_out_Lpwr(i_part,i_elect,i_f,i_t) = model_out_Lpwr_trl(2);
                
                clear mdn_pwr model_out_Hpwr_trl model_out_Lpwr_trl
            end
            clear i_f
        end
        clear i_t tmp_pwr i_elect 
    end
    clear ii i_elect
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
%% Save Data
% this is large file so it will take some time to save
save([saveLocation 'mdn_pwrR_split_v4.mat'],...
    'errdeg_Hpwr_trl','errdeg_Lpwr_trl','g_out_Hpwr','g_out_Lpwr',...
    'sd_out_Hpwr','sd_out_Lpwr','model_out','times','-v7.3')

% /////////////////////////////////////////////////////////////////////////










