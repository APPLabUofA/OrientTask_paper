% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Need to have loaded the single trial data using LoadProcData_OrientTask.m
% Should end up with all_ersp cell variable
% Code below is for regular trials (not catch trials) that has been aligned
% to target onset

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat'); %setting for target-aligned trials except catch trials

%% Location to save power data
saveLocation = [exp.dataLocation '\ProcessData\']; % set save directory of data set

%% Location to save figures
saveFig = [pwd '\Figures\ERS\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

%% Load saved behavioral data
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% Raw ERS values (log scaled)
% /////////////////////////////////////////////////////////////////////////
% #########################################################################

% --For data with targets--
all_erspN = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % --
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        tmp_ersp = abs(all_ersp{i_part,i_elect}).^2;
        for i_trial = 1:size(tmp_ersp,3)
            all_erspN{i_part,i_elect}.trials(:,:,i_trial) = log10(tmp_ersp(:,:,i_trial)); %dB converted
        end
        clear i_trial
    end
    clear ii i_elect tmp_ersp
end
clear i_part



% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%%                      Standardize Power
% /////////////////////////////////////////////////////////////////////////
% #########################################################################

%finds the times you want from the timess variable
timewin = [-700 567];
timephi = find(times>=timewin(1) & times<=timewin(2));

all_ersp_Z = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
% Change power to z-score values per person
for i_part = 1:length(exp.participants)
    % Get power across trials
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_erspN{i_part,i_elect}.trials; %get single subject's baseline corrected power
%         all_ersp_Z{i_part,i_elect}.trials = normalize(part_ersp,3,'zscore','robust');
        
        baseline_power = part_ersp(:,timephi,:);
        baselineZ = (part_ersp-repmat(mean(baseline_power,2),1,size(part_ersp,2))) ./ repmat(std(baseline_power,[],2),1,size(part_ersp,2));
        all_ersp_Z{i_part,i_elect}.trials = baselineZ;
%         all_ersp_Z{i_part,i_elect}.trials =...
%             (part_ersp - mean(part_ersp(:,timephi,:),2)) ./ std(part_ersp(:,timephi,:),[],2);
        clear part_ersp i_elect baselineZ baseline_power
    end
    clear ii
end
clear i_part timephi


% /////////////////////////////////////////////////////////////////////////
%% Save Standardized Data
chanlocs = EEG.chanlocs; %going to want to save electrode locations
% this is large file so it will take some time to save
save([saveLocation 'all_ersp_Z_v4.mat'],'all_ersp_z','chanlocs','timewin',...
    '-v7.3')


% /////////////////////////////////////////////////////////////////////////
%% OR Load Standardized Data If Exists
all_ersp_Z = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'all_ersp_z'));  %gets loaded as a struct
all_ersp_Z = all_ersp_Z{1};
chanlocs = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};
timewin = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'timewin')); 
timewin = timewin{1};

% /////////////////////////////////////////////////////////////////////////




% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% Plots in Supporting Info Figures
% /////////////////////////////////////////////////////////////////////////
% #########################################################################



% /////////////////////////////////////////////////////////////////////////
%% # ERS: Power by Model SD #
% /////////////////////////////////////////////////////////////////////////

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Load previously created data (if it has been created)
load([saveLocation 'pwr_AvG_v4.mat']);
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


%% Create ERS by errors
x_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
n_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
x_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
n_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
errlims = cell(1,length(exp.participants));    %pre-allocate
for i_part = 1:length(exp.participants)
    
    % Get upper and lower limits based on model fit
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    % Get errors values
    x_errdeg_m{i_part} = resp_errdeg{i_part}(resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)); %small errors
    n_errdeg_m{i_part} = resp_errdeg{i_part}([find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))]);

    % Calculate power
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp_Z{i_part,i_elect}.trials; %get single subject's baseline corrected power
        
        % Get trials with small errors
        x_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[...
            find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] ),3));

        % Get trials with large errors
        n_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[...
            find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] ),3));
        
        clear part_ersp i_elect
    end
end
clear ii i_part

% /////////////////////////////////////////////////////////////////////////
% Save data if not saved yet
save([saveLocation 'pwr_AvG_v4.mat'],'errlims','n_errdeg_m','n_pwr',...
    'x_errdeg_m','x_pwr','timewin');

% /////////////////////////////////////////////////////////////////////////



% /////////////////////////////////////////////////////////////////////////
% #########################################################################
%% Gets a count of trials
err_trl_count(:,1) = cellfun(@numel,x_errdeg_m); %small errors
err_trl_count(:,2) = cellfun(@numel,n_errdeg_m); %large errors
% err_trl_count(:,3) = cell2mat({ALLEEG(1:end).trials}); %total trial count
% #########################################################################
% /////////////////////////////////////////////////////////////////////////
% #########################################################################


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% Plot spectogram across subjects &&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Raw ERS plots
cmap = jet; %create colormap colors
savename = 'SpecPlotZ_';
for ii = 1:length(exp.singletrialselecs)
% for ii = 1:5 %central electrodes only
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    
    CLim = [-0.4 0.4]; %set power scale of plot
    
    % Plot Small Errors
    figure('Position', [1 1 1685 405]); colormap(cmap) %open a new figure
    subplot(1,2,1)
    imagesc(times,freqs,plot_ers_x,CLim);
    title(['Accurate: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Normalized Power');
    
    % Plot Large Errors
    subplot(1,2,2)
    imagesc(times,freqs,plot_ers_n,CLim);
    title(['Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Normalized Power');
    
    savefig([saveFig savename exp.singtrlelec_name{ii}])
   
    clear plot_ers_x plot_ers_n CLim t
end
clear ii i_elect cmap savename

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Difference ERS plot
cmap = redblue(256); %create colormap colors
savename = 'SpecPlot_DifZ_';
for ii = 1:length(exp.singletrialselecs)
% for ii = 1:5
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    
    CLim = [-0.3 0.3]; %set power scale of plot
    
    % Plot Accurate-Guesses
    figure; colormap(cmap) %open a new figure
    imagesc(times,freqs,plot_ers_x-plot_ers_n,CLim);
    title(['Accurate-Guesses: ' exp.singtrlelec_name{ii}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca); 
    t.Ticks = [-0.3:0.1:0.3]; %make sure colorbar contains ticks
    set(get(t,'ylabel'),'String', 'Normalized Power Difference');
    
    savefig([saveFig savename exp.singtrlelec_name{ii}])
    
    clear plot_ers_x plot_ers_n CLim t
end
clear ii i_elect cmap savename


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% Grand average ERS plot

plot_ers_x = NaN(length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
plot_ers_n = NaN(length(exp.singletrialselecs),length(freqs),length(times)); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x(i_elect,:,:) = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n(i_elect,:,:) = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
end
clear ii i_elect

%Grand Average
plot_x_avg = squeeze(nanmean(plot_ers_x,1)); %small errors
plot_n_avg = squeeze(nanmean(plot_ers_n,1)); %large errors
clear plot_ers_x plot_ers_n

CLim = [-0.4 0.4]; %set power scale of plot
cmap = jet; %create colormap colors

% Open new figure
figure('Position', [1 1 1685 405]); colormap(cmap)

% Plot Small Errors
subplot(1,2,1)
imagesc(times,freqs,plot_x_avg,CLim);
title('Accurate: Grand Avg'); set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Normalized Power');
    
% Plot Large Errors
subplot(1,2,2)
imagesc(times,freqs,plot_n_avg,CLim);
title('Guesses: Grand Ave'); set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Normalized Power');

savefig([saveFig 'SpecPlotZ_GrandAvg'])
   
clear CLim t cmap
    

% .........................................................................
% Plot Accurate-Guesses
CLim = [-0.3 0.3]; %set power scale of plot
cmap = redblue(256); %create colormap colors
figure; colormap(cmap) %open a new figure
imagesc(times,freqs,plot_x_avg-plot_n_avg,CLim);
title('Accurate-Guesses: Grand Avg'); set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
% xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca); 
t.Ticks = [-0.3:0.1:0.3]; %make sure colorbar contains ticks
set(get(t,'ylabel'),'String', 'Normalized Power Difference');

savefig([saveFig 'SpecPlot_DifZ_GrandAvg'])

clear plot_x_avg plot_n_avg CLim t cmap

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     ''''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% List electrodes to get ERP topograph plots (need all of them) 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

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
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

% Get mean power at frequency band for each electrode
pwr_top = NaN([length(times),length(elect_erp),2]); %pre-allocate
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
                                                 %{electrode}(part,freq,time)
    pwr_top(:,i_elect,1) = squeeze(mean(mean(x_pwr{1,i_elect}(:,freqlim,:),2),1)); %small errors
    pwr_top(:,i_elect,2) = squeeze(mean(mean(n_pwr{1,i_elect}(:,freqlim,:),2),1)); %large errors
end
clear ii i_elect

CLim = [-3 3]; %set power scale of plot
CLim2 = [-0.3 0.3]; %set power scale of plot difference
colormap('jet')
for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %finds the times you want from the times variable
    time_window = find(times>= itWin(1),1):find(times>= itWin(2),1)-1;
    
    % ---Plot trial type power maps---
    figure('Color',[1 1 1],'Position',[1 1 941 349]);
%     figure('Color',[1 1 1]);
    set(gca,'Color',[1 1 1]);
    
    temp = mean(pwr_top(time_window,:,1),1)'; %get power at time window
    temp(1) = NaN; %not M2 electrode
    
    subplot(1,2,1)
    topoplot(temp,chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})
    title('Accurate');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Standardized Power (dB)');
    clear temp
    hold on
    
    temp = mean(pwr_top(time_window,:,2),1)'; %get power at time window
    temp(1) = NaN; %not M2 electrode
    
    subplot(1,2,2)
    topoplot(temp,chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})
    title('Guesses');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Standardized Power');
    clear temp
    
    supertitle([num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms']);
    
    savefig([saveFig 'Topo_' num2str(freqband(1)) '-' num2str(freqband(2)) '_' num2str(itWin(1)) 'to' num2str(itWin(2))])
    
    hold off
    
    % ---Plot difference power maps---
    figure('Color',[1 1 1]);
    set(gca,'Color',[1 1 1]);
    
    temp = (mean(pwr_top(time_window,:,1),1)-mean(pwr_top(time_window,:,2),1))'; %get difference in power at time window
    temp(1) = NaN; %not M2 electrode
    topoplot(temp,chanlocs,'whitebk','on','plotrad',0.6,'maplimits',CLim2,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})
    title(['Accurate-Guesses: ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms']);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Standardized Power');
    clear temp
    
    savefig([saveFig 'TopoDiff_' num2str(freqband(1)) '-' num2str(freqband(2)) '_' num2str(itWin(1)) 'to' num2str(itWin(2))])

    
    clear itWin time_window temp
end
clear tw_i t

clear freqlim freqband CLim pwr_top CLim2



clear tWin






