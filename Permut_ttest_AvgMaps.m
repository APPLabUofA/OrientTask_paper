% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 
currentfolder = pwd; %to return to current folder after loaded data

%% Location of saved permutation results
%  Folder location should contain the repeated permutation test results
%below is path for data created by Permut_Paired_ttest_Pwr.m
saveLocation_perm = [saveLocation 'permut_data\paired_ttest_pwr\']; 

%% Location to save figures
%usually save in the same folder as the permuted results 
saveFig = [saveLocation_perm 'Figures\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load saved mat files from permutation testing code

%*Change to folder with saved permutation test files before running code below*
cd(saveLocation_perm)


zmaplist = cellstr(ls('*.mat')); %get list of saved data in folder

% load data and put into variables
combo_threshold = cell(1,length(zmaplist)); %pre-allocate
combo_zmap = cell(1,length(zmaplist)); %pre-allocate
for jj = 1:length(zmaplist)
    load(zmaplist{jj})
    combo_threshold{1,jj} = threshold_out;
    combo_zmap{1,jj} = zmap_out;
    clear threshold_out zmap_out
end
clear jj zmaplist


for ii = 1:length(combo_threshold)
    for ip = 1:length(exp.singletrialselecs) 
        i_elect = exp.electrode(ip); %get electrode number
        
        c_thresh(ii,i_elect,1:2) = combo_threshold{1,ii}{i_elect};
        c_zmap(ii,i_elect,:,:) = combo_zmap{1,ii}{i_elect};
    end
    clear i_elect
end
clear ii ip


% Get mean of maps and threshold
cout_thresh = squeeze(mean(c_thresh,1));
cout_zmap = squeeze(mean(c_zmap,1));
clear c_thresh c_zmap

% Return to main folder
cd(currentfolder)


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Plot averaged permutation results

cmap = redblue(256); %create colormap colors
CLim = [-6 6]; %color axis limits

figname = 'Permut_AvgMap_pwr_'; %name of saved figure

for ii = 1:length(exp.singletrialselecs)
% for ii = 1:5

    i_elect = exp.electrode(ii); %get electrode number
    
    % Get threshold values
    lower_threshold = cout_thresh(i_elect,1);
    upper_threshold = cout_thresh(i_elect,2);
    
    % Plot the uncorrected z-values
    figure; colormap(cmap)
    contourf(times,freqs,squeeze(cout_zmap(i_elect,:,:)),40,'linecolor','none')
    
    % Apply corrected threshold contour map
    zmapthresh = squeeze(cout_zmap(i_elect,:,:));
    zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=false;
    zmapthresh=logical(zmapthresh);
    hold on
    contour(times,freqs,zmapthresh,1,'linecolor','k') %plots black lines

    axis square
    set(gca,'clim',CLim)
    title(['Thresholded Z map: ' exp.singtrlelec_name{ii}],'FontSize',14);  
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5) %vertical line for response screen onset
    ylim([2 40]); yticks(5:5:40)
%     xlim([-700 800]); xticks(-600:200:800)
    xlim([-200 800]); xticks(-200:100:800) %to match ERPs
    ylabel('Frequency (Hz)');
    xlabel('Time (ms)');
    t = colorbar('peer',gca);
    t.Label.String = 'Z-Score';
    
    savefig([saveFig figname exp.elec_names{ii}])
    
    clear zmapthresh lower_threshold upper_threshold
    
end
clear ii i_elect t cmap figname 



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Plot Grand Average Permutation Results

cmap = redblue(256); %create colormap colors
CLim = [-6 6]; %color axis limits

figname = 'Permut_AvgMap_pwr_'; %name of saved figure

% Apply Corrected threshold contour map
lower_threshold = mean(cout_thresh(:,1),1);
upper_threshold = mean(cout_thresh(:,2),1);

figure; colormap(cmap)
contourf(times,freqs,squeeze(mean(cout_zmap(:,:,:),1)),40,'linecolor','none')
zmapthresh = squeeze(mean(cout_zmap(:,:,:),1));
zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=false;
zmapthresh=logical(zmapthresh);
hold on
contour(times,freqs,zmapthresh,1,'linecolor','k')

axis square
set(gca,'clim',CLim)
title('Thresholded Z map: Grand Average','FontSize',14);  
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5) %vertical line for response screen onset
ylim([3 40]); yticks(5:5:40)
xlim([-200 800]); xticks(-200:100:800) %axis limit to match ERPs
% xlim([-700 800]); xticks(-600:200:800)
ylabel('Frequency (Hz)');
xlabel('Time (ms)');
t = colorbar('peer',gca);
t.Label.String = 'Z-Score';

savefig([saveFig figname 'GrandAvg'])

clear zmapthresh lower_threshold upper_threshold
    

clear ii i_elect t cmap figname 


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%don't clear if doing topographys
clear combo_threshold combo_zmap cout_thresh cout_zmap

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     ''''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% These plots were not included in the paper, but can be useful to look at

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% Need structure with electrode locations
chanlocs = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};

% Set the range of time to consider
tWin{1} = [0 80];
tWin{2} = [80 140]; %P1
tWin{3} = [140 200]; %N1
tWin{4} = [200 255]; %P2
tWin{5} = [255 360]; %N2
tWin{6} = [360 500]; %P3

%finds the frequencies you want
freqband = [30 40]; %gamma
% freqband = [23 29]; %beta2
% freqband = [15 22]; %beta1
% freqband = [10 14]; %high alpha
% freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

% Get mean power at frequency band for each electrode
pwr_top = NaN([length(times),length(elect_erp),1]); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    pwr_top(:,i_elect,1) = squeeze(mean(cout_zmap(i_elect,freqlim,:),2));
end
clear ii i_elect

cmap = redblue(256); %create colormap colors
CLim = [-6 6]; %color axis limits
figname = 'Perm_AvgPwrTop_';

for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %finds the times you want from the times variable
    time_window = find(times>= itWin(1),1):find(times>= itWin(2),1)-1;

    temp = mean(pwr_top(time_window,:,1),1)'; %get power at time window
    temp(1) = NaN; %not M2 electrode
    
    %make figure with white background
    figure('Color',[1 1 1]); set(gca,'Color',[1 1 1]);
    
    topoplot(temp,chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})
    colormap(cmap)
    title(['Accurate v Guess: ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms']);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Standardized Power');

    savefig([saveFig figname num2str(freqband(1)) '-' num2str(freqband(2)) '_' num2str(itWin(1)) 'to' num2str(itWin(2))])
    
    clear itWin time_window temp
end
clear tw_i t cmap

clear freqlim freqband CLim pwr_top


clear tWin figname 
close all

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

clear combo_threshold combo_zmap cout_thresh cout_zmap elect_erp el_erp_names

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

