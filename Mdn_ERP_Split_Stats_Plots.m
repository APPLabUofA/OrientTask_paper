% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% This code uses data previously processed by Mdn_ERP_Split.m
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
% Load data from Mdn_ERP_Split.m
load([saveLocation 'mdn_ERP_split_v4.mat'])

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% /////////////////////////////////////////////////////////////////////////
%% Location to save figures
saveFigLoc = [pwd '\Figures\mdn_ERP_split\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFigLoc)
    mkdir(saveFigLoc);
end


% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% :::::::::::::::::::::  Permutation Test  :::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 
% **Need code from: https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox
% ------------------------------------------------------------------------- 

%% Guess Rate

alpha = 0.05; %Set alpha level
nperm = 1e5; %Number of permutations

% Re-order dimensions
statA(:,:,1) = g_aboveP1';
statA(:,:,2) = g_aboveN1';
statA(:,:,3) = g_aboveP2';
statA(:,:,4) = g_aboveN2';
statA(:,:,5) = g_aboveP3';
A = permute(statA,[1 3 2]);

statB(:,:,1) = g_belowP1';
statB(:,:,2) = g_belowN1';
statB(:,:,3) = g_belowP2';
statB(:,:,4) = g_belowN2';
statB(:,:,5) = g_belowP3';
B = permute(statB,[1 3 2]);

% One sample/repeated-measures permutation test
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(A-B),nperm,alpha);

% Average across subjects
mean_g_all(:,:,1) = squeeze(mean(statA,2));
std_g_all(:,:,1) = squeeze(std(statA,[],2));
mean_g_all(:,:,2) = squeeze(mean(statB,2));
std_g_all(:,:,2) = squeeze(std(statB,[],2));

clear statA statB pval t_orig tmx_ptile nperm alpha A B


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

%% SD Parameter

alpha = 0.05; %Set alpha level
nperm = 1e5; %Number of permutations

% Re-order dimensions
statA(:,:,1) = sd_aboveP1';
statA(:,:,2) = sd_aboveN1';
statA(:,:,3) = sd_aboveP2';
statA(:,:,4) = sd_aboveN2';
statA(:,:,5) = sd_aboveP3';
A = permute(statA,[1 3 2]);

statB(:,:,1) = sd_belowP1';
statB(:,:,2) = sd_belowN1';
statB(:,:,3) = sd_belowP2';
statB(:,:,4) = sd_belowN2';
statB(:,:,5) = sd_belowP3';
B = permute(statB,[1 3 2]);

% One sample/repeated-measures permutation test
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(A-B),nperm,alpha);


% Average across subjects
mean_sd_all(:,:,1) = squeeze(mean(statA,2));
std_sd_all(:,:,1) = squeeze(std(statA,[],2));
mean_sd_all(:,:,2) = squeeze(mean(statB,2));
std_sd_all(:,:,2) = squeeze(std(statB,[],2));


clear statA statB pval t_orig tmx_ptile nperm alpha A B 



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


%% Guess Rate

% Difference Between Conditions
mean_g_all(:,:,3) = mean_g_all(:,:,1) - mean_g_all(:,:,2);

nconds = 3; %plotting the different conditions
conds = {'Above Median';'Below Median';'Above - Below'}; %labels for plots
ERP_labels = {'P1';'N1';'P2';'N2';'P3'}; %supertitle

for ierp = 1:size(mean_g_all,2) %loop through ERPs
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);
    tiledlayout(1,3)

    for i_cond = 1:nconds
        
        temp = NaN(length(elect_erp),1);
        temp(2:32,1) = squeeze(mean_g_all(:,ierp,i_cond));
        
        if i_cond == 3
            CLims = [-0.1 0.1];
%             subplot(1,3,i_cond); 
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar labels
%             colorbar('Ticks',[-0.09,-0.06,-0.03,0,0.03],...
%                 'TickLabels',{'-0.09','-0.06','-0.03','0.00','0.03'})
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Guess Rate (g) Difference');
            clear CLims
        else
            CLims = [0.12 0.22];
%             subplot(1,3,i_cond);
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
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
    supertitle(['Trials Split By ' ERP_labels{ierp} ' Amplitude'],...
        'FontSize',10.5)
    
    savefig([saveFigLoc 'g_topo_' ERP_labels{ierp} '.fig']); %save as matlab figure
    saveas(gcf,[saveFigLoc 'g_topo_' ERP_labels{ierp}],'svg'); %save for adobe illustrator

    clear i_cond t
end
clear ierp nconds conds ERP_labels saveFigLoc


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% SD

% Difference Between Conditions
mean_sd_all(:,:,3) = mean_sd_all(:,:,1) - mean_sd_all(:,:,2);

nconds = 3; %plotting the different conditions
conds = {'Above Median';'Below Median';'Above - Below'}; %labels for plots
ERP_labels = {'P1';'N1';'P2';'N2';'P3'}; %supertitle

for ierp = 1:size(mean_sd_all,2) %loop through ERPs
    
    figure('Color',[1 1 1],'Position',[1 1 941 349]);
    tiledlayout(1,3)

    for i_cond = 1:nconds
        
        temp = NaN(length(elect_erp),1);
        temp(2:32,1) = squeeze(mean_sd_all(:,ierp,i_cond));
        
        if i_cond == 3
            CLims = [-1.5 1.5];
%             subplot(1,3,i_cond); 
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
                'plotchans',elect_erp,'emarker',{'.','k',11,1})
            % Color bar labels
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Standard Deviation (SD) Difference');
            clear CLims
        else
            CLims = [9 12];
%             subplot(1,3,i_cond);
            nexttile 
            set(gca,'Color',[1 1 1]);
            topoplot(temp,EEG.chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLims,...
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
    supertitle(['Trials Split By ' ERP_labels{ierp} ' Amplitude'],...
        'FontSize',10.5)
    
    savefig([saveFigLoc 'sd_topo_' ERP_labels{ierp} '.fig']); %save as matlab figure
    saveas(gcf,[saveFigLoc 'sd_topo_' ERP_labels{ierp}],'svg'); %save for adobe illustrator

    clear i_cond t
end
clear ierp nconds conds ERP_labels saveFigLoc



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% clears variables that end/begin with...
clear -regexp \<permtest_ \<mean_ \<std_ \<g_ \<sd_ \<h_ \<errdeg_ \<crit_...
    \<adj_
clear timelim nperm alpha time_points elect_erp el_erp_names h time_win

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ------------------------------------------------------------------------- 



















