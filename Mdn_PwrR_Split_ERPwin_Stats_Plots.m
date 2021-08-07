% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% This code uses data previously processed by Mdn_PwrR_Split.m
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
% Load data saved by Mdn_PwrR_Split.m
load('mdn_pwrR_split_v4.mat')

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% /////////////////////////////////////////////////////////////////////////
%% Location to save figures
saveFigLoc = [pwd '\Figures\mdn_pwr_split_v4\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFigLoc)
    mkdir(saveFigLoc);
end



% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%%            Guess Rate in Time and Frequency Windows
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

clear timewin freqlim freqband

%finds the frequencies you want
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

g_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
g_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2));
                                                        % (participant x electrode x freq x time)
            g_Hpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(g_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
            g_Lpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(g_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            clear timelim
        end
    end
    clear ii i_elect
end
clear i_part


% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% :::::::::::::::::::::  Permutation Test  :::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 
% **Need code from: https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox
% ------------------------------------------------------------------------- 

alpha = 0.05/2; %Set alpha level (bonferroni correction for freq bands)
nperm = 1e5; %Number of permutations

%data - 3D matrix of ERPs from Group A (Channel x Time x Participant)
statA = permute(g_Hpwr_stat,[2 3 1]); % re-order dimensions
statB = permute(g_Lpwr_stat,[2 3 1]); % re-order dimensions

[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(statA-statB),nperm,alpha);


clear statA statB t_orig tmx_ptile pval g_Hpwr_stat g_Lpwr_stat timewin...
    freqlim freqband




% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%%               SD in Time and Frequency Windows
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

clear timewin freqlim freqband

%finds the frequencies you want
freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

sd_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
sd_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2));
                                                        % (participant x electrode x freq x time)
            sd_Hpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(sd_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
            sd_Lpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(sd_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            clear timelim
        end
    end
    clear ii i_elect
end
clear i_part



% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% :::::::::::::::::::::  Permutation Test  :::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 
% **Need code from: https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox
% ------------------------------------------------------------------------- 

alpha = 0.05/2; %Set alpha level
nperm = 1e5; %Number of permutations

%data - 3D matrix of ERPs from Group A (Channel x Time x Participant)
statA = permute(sd_Hpwr_stat,[2 3 1]); % re-order dimensions
statB = permute(sd_Lpwr_stat,[2 3 1]); % re-order dimensions

[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(statA-statB),nperm,alpha);


clear statA statB t_orig tmx_ptile pval sd_Hpwr_stat sd_Lpwr_stat timewin...
    freqlim freqband






% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''    Guess Rate: Topographys     '''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% /////////////////////////////////////////////////////////////////////////

clear timewin freqlim freqband

%finds the frequencies you want
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

g_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
g_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2));
                                                        % (participant x electrode x freq x time)
            g_Hpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(g_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
            g_Lpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(g_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            clear timelim
        end
    end
    clear ii i_elect
end
clear i_part

% /////////////////////////////////////////////////////////////////////////

% Mean across subjects
mean_g_all(:,:,1) = squeeze(mean(g_Hpwr_stat,1));
mean_g_all(:,:,2) = squeeze(mean(g_Lpwr_stat,1));
mean_g_all(:,:,3) = mean_g_all(:,:,1) - mean_g_all(:,:,2); % Difference
clear g_Hpwr_stat g_Lpwr_stat

% /////////////////////////////////////////////////////////////////////////


nconds = 3; %plotting the different conditions
conds = {'High Log Power';'Low Log Power';'Difference (H-L)'}; %labels for plots
ERP_labels = {'80-140 ms (P1)';'140-200 ms (N2)';'200-255 ms (P2)';...
    '255-360 ms (N2)';'360-500 ms (P3)'}; %supertitle


CLims1 = [0.12 0.22];
CLims3 = [-0.1 0.1];
% CLims1 = [0.14 0.21];
% CLims3 = [-0.04 0.02];

for ierp = 1:length(timewin) %loop through ERPs
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);
    tiledlayout(1,3)

    for i_cond = 1:nconds
        
        temp = NaN(length(elect_erp),1);
        temp(2:32,1) = squeeze(mean_g_all(:,ierp,i_cond));
        
        if i_cond == 3
%             subplot(1,3,i_cond); 
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims3,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar labels
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Guess Rate (g) Difference');
            clear CLims
        else
%             subplot(1,3,i_cond);
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims1,...
            'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar label
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Guess Rate (g)');
            clear CLims
        end
        
        title(conds{i_cond});
        clear temp t
    end
    
    % Overall subplot title
    supertitle([num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz Frequency: '...
        ERP_labels{ierp}],'FontSize',10.5)
    
    % Save Figure
    savename = ['g_topo_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_'...
        num2str(timewin{ierp}(1)) '-' num2str(timewin{ierp}(2)) 'ms']; %name to save figures
    savefig([saveFigLoc savename '.fig']); %save as matlab figure
    saveas(gcf,[saveFigLoc savename],'svg'); %save for adobe illustrator

    clear i_cond t savename
end
clear ierp nconds conds ERP_labels saveFigLoc freqband timewin mean_g_all...
    CLims1 CLims3



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''    SD: Topographys     '''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% /////////////////////////////////////////////////////////////////////////

clear timewin freqlim freqband

%finds the frequencies you want
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

sd_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
sd_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2));
                                                        % (participant x electrode x freq x time)
            sd_Hpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(sd_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
            sd_Lpwr_stat(i_part,ii,iwin) = squeeze(mean(mean(sd_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            clear timelim
        end
    end
    clear ii i_elect
end
clear i_part

% /////////////////////////////////////////////////////////////////////////

% Mean across subjects
mean_sd_all(:,:,1) = squeeze(mean(sd_Hpwr_stat,1));
mean_sd_all(:,:,2) = squeeze(mean(sd_Lpwr_stat,1));
mean_sd_all(:,:,3) = mean_sd_all(:,:,1) - mean_sd_all(:,:,2); % Difference
clear sd_Hpwr_stat sd_Lpwr_stat

% /////////////////////////////////////////////////////////////////////////


nconds = 3; %plotting the different conditions
conds = {'High Log Power';'Low Log Power';'Difference (H-L)'}; %labels for plots
ERP_labels = {'80-140 ms (P1)';'140-200 ms (N2)';'200-255 ms (P2)';...
    '255-360 ms (N2)';'360-500 ms (P3)'}; %supertitle

CLims1 = [9 12];
CLims3 = [-1.5 1.5];

for ierp = 1:length(timewin) %loop through ERPs
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);
    tiledlayout(1,3)

    for i_cond = 1:nconds
        
        temp = NaN(length(elect_erp),1);
        temp(2:32,1) = squeeze(mean_sd_all(:,ierp,i_cond));
        
        if i_cond == 3
%             subplot(1,3,i_cond); 
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims3,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar labels
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Standard Deviation (SD) Difference');
            clear CLims
        else
%             subplot(1,3,i_cond);
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims1,...
            'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar label
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Standard Deviation (SD)');
            clear CLims
        end
        
        title(conds{i_cond});
        clear temp t
    end
    
    % Overall subplot title
    supertitle([num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz Frequency: '...
        ERP_labels{ierp}],'FontSize',10.5)
    
    % Save Figure
    savename = ['sd_topo_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_'...
        num2str(timewin{ierp}(1)) '-' num2str(timewin{ierp}(2)) 'ms']; %name to save figures
    savefig([saveFigLoc savename '.fig']); %save as matlab figure
    saveas(gcf,[saveFigLoc savename],'svg'); %save for adobe illustrator

    clear i_cond t savename
end
clear ierp nconds conds ERP_labels saveFigLoc freqband timewin mean_sd_all...
    CLims1 CLims3


% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////















