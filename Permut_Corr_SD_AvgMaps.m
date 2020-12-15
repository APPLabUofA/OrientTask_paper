
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
%below is path for data created by Permut_Corr_SD_Pwr.m
saveLocation_perm = [saveLocation 'permut_data\Corr_SD_Pwr\']; 

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
combo_realcorrs = cell(1,length(zmaplist)); %pre-allocate
combo_realcorrs_plot = cell(1,length(zmaplist)); %pre-allocate
for jj = 1:length(zmaplist)
    load(zmaplist{jj})
    combo_threshold{1,jj} = threshold_out;
    combo_zmap{1,jj} = zmap_out;
    combo_realcorrs{1,jj} = realcorrs_out;
    combo_realcorrs_plot{1,jj} = realcorrs_plot_out;
    clear threshold_out zmap_out realcorrs_out realcorrs_plot_out
end
clear jj zmaplist

for ii = 1:length(combo_threshold)
    for ip = 1:length(exp.singletrialselecs) 
%     for ip = 1 %grand mean    
        i_elect = exp.electrode(ip); %get electrode number
        
        c_thresh(ii,i_elect,1:2) = combo_threshold{1,ii}{ip};
        c_zmap(ii,i_elect,:,:) = combo_zmap{1,ii}{ip};
        c_realcorrs(ii,i_elect,:,:) = combo_realcorrs{1,ii}{ip};
        c_realcorrs_plot(ii,i_elect,:,:) = combo_realcorrs_plot{1,ii}{ip};
        
    end
    clear i_elect
end
clear ii ip

% Get mean of maps and threshold
cout_thresh = squeeze(mean(c_thresh,1));
cout_zmap = squeeze(mean(c_zmap,1));
cout_realcorrs = squeeze(mean(c_realcorrs,1));
cout_realcorrs_plot = squeeze(mean(c_realcorrs_plot,1));
clear c_thresh c_zmap c_realcorrs c_realcorrs_plot


% Return to main folder
cd(currentfolder)

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Plot averaged permutation results

cmap = redblue(256); %create colormap colors
CLim = [-1 1]; %color axis limits

figname = 'Avg_Corr_SD_'; %name of saved figure

for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.electrode(ii); %get electrode number

    % apply Corrected threshold contour map
    lower_threshold = cout_thresh(i_elect,1);
    upper_threshold = cout_thresh(i_elect,2);
    
    % open new figure
    figure; colormap(cmap)
    contourf(times,freqs,squeeze(cout_realcorrs_plot(i_elect,:,:)),40,'linecolor','none')
    
    zmapthresh = squeeze(cout_zmap(i_elect,:,:));
    realcorrs = squeeze(cout_realcorrs(i_elect,:,:));  
    zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=false;
    
    zmapthresh=logical(zmapthresh);
    hold on
    contour(times,freqs,zmapthresh,1,'linecolor','k')

    axis square
    set(gca,'clim',CLim)
    title(['Spearman Rho map: ' exp.singtrlelec_name{ii}],'FontSize',14);  
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5) %vertical line for response screen onset
%     ylim([3 40]); yticks(5:5:40)
%     xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
    xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Frequency (Hz)');
    xlabel('Time (ms)');
    t = colorbar('peer',gca);
    t.Label.String = 'Spearman"s Rho';
    
    savefig([saveFig figname exp.elec_names{ii}])
    
    clear zmapthresh lower_threshold upper_threshold

 
end
clear ii i_elect t cmap figname CLim 



clear combo_threshold combo_zmap cout_thresh cout_zmap combo_realcorrs...
    combo_realcorrs_plot cout_realcorrs cout_realcorrs_plot saveFig



% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Scatter Plots ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% Re-create scatter plots in manuscript

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%% Load previously processed target-aligned epoch data
%data created in Power_ZScore_ERS_Plots.m
all_ersp_Z = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'all_ersp_z'));  %gets loaded as a struct
all_ersp_Z = all_ersp_Z{1};

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Extract standard mixed model parameter values
% loaded from saved behavioral data or can be created in BEH_ModelFits_OrientTask.m
sd_out_cat = NaN(length(exp.participants),1);
% g_out_cat = NaN(length(exp.participants),1);
for i_part = 1:length(exp.participants)
%     g_out_cat(i_part,1) = model_out{1,i_part}(1);
    sd_out_cat(i_part,1) = model_out{1,i_part}(2);
end
clear i_part

tmp_out = sd_out_cat;

% tmp_out(8,:) = []; %remove subject 8 (outlier in frontal electrodes)


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Reorganize time-frequency data
% obs_pwr{electrode}(freqs x time points x subjects)
obs_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    for i_part = 1:length(exp.participants)
        obs_pwr{1,i_elect}(:,:,i_part) = squeeze(mean(all_ersp_Z{i_part,i_elect}.trials,3)); %get single subject's z power
    end
    clear i_part
end
clear ii i_elect

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% -------------------------------------------------------------------------
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Select locations and times of interest

clear timewin freqlim freqband

%finds the frequencies you want
% freqband = [14 17]; %PO4
% freqband = [14 17]; %O2
% freqband = [19 22]; %O2-2
% freqband = [15 21]; %O2-3
% freqband = [13 18]; %PO4
freqband = [2.8 3.2]; %F4 & F3
% freqband = [14 15]; %P6
% freqband = [15 16]; %Oz
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

%finds the times you want from the timess variable
% timewin = [190 270];%PO4---
% timewin = [215 280];%O2
% timewin = [250 320];%O2-2
% timewin = [380 440];%O2-3
% timewin = [150 300];%PO4
timewin = [535 590];%F4 (55 ms)
% timewin = [555 595];%F3 (40 ms)
% timewin = [230 260];%P6
% timewin = [230 250];%Oz
timelim = find(times>=timewin(1) & times<=timewin(2));

% Select electrode
% ii = 9; %PO4
% ii = 7; %O2
% ii = 1; %Oz
% ii = 13; %P6
ii = 29; %F4
% ii = 28; %F3
i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
disp(exp.elec_names{ii}) %make sure selected correct electrode

plotpwr = squeeze(mean(mean(obs_pwr{1,i_elect}(freqlim,timelim,:),2),1));
% plotpwr(8,:) = []; %remove subject 8


%% Bootstrap the correlations
scorr = @(a,b)(corr(a,b,'type','Spearman')); %to use spearman's rho correlation
n_samples = 1000;  %This is the number of random samples for bootstrapping, make it smaller if it takes too long, bigger for nicer pics

[cstat cids] = bootstrp(n_samples,scorr,plotpwr,tmp_out); 
% [cstat cids] = bootstrp(n_samples,@corr,plotpwr,tmp_out);
boot_mean = mean(cstat);
boot_std_err = std(cstat);
p_corr_boot = length(find(cstat<0))/n_samples;
%Calculate CI with bootstrapping
ci = bootci(n_samples,scorr,plotpwr,tmp_out);
% ci = bootci(n_samples,@corr,plotpwr,tmp_out);


%% Plot each correlation line 
% Make this more than 10 but less than 100 or it will take forever
n_lines = 60;
figure
for i_sample=1:n_lines
    scatter(plotpwr(cids(:,i_sample),1),tmp_out(cids(:,i_sample),1));
    hold on
    lsline;
end
scatter(plotpwr,tmp_out,'MarkerFaceColor',[0 .7 .7],'MarkerEdgeColor','k')
ylabel('SD (\sigma)','Interpreter','tex'); xlabel('Standardized Power (Z)');
ylim([4 20])
title([exp.elec_names{ii}...
    ': Mean rho = ' num2str(round(boot_mean,2)) ', 95% CI = [' num2str(round(ci(1),2))...
    ' '  num2str(round(ci(2),2)) '], p = ' num2str(round(p_corr_boot,2))]);

% xlim([-0.1 0.4]); 
% xticks(-0.1:0.1:0.4)

hold off


clear timewin freqlim freqband plotpwr freqband timelim ii i_elect n_lines cstat cids...
    boot_std_err boot_mean ci p_corr_boot











