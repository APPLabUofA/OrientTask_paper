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


%% Load raw power data
% Data name is for dataset provided. This name will be different if data is 
% processed again
load([pwd 'all_ersp_v4.mat'])

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%%                     Raw ERS Values (log scaled)
% /////////////////////////////////////////////////////////////////////////
% #########################################################################
% --For data with targets--
all_erspR = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % --
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        tmp_ersp = abs(all_ersp{i_part,i_elect}).^2;
        for i_trial = 1:size(tmp_ersp,3)
            all_erspR{i_part,i_elect}.trials(:,:,i_trial) = log10(tmp_ersp(:,:,i_trial)); %dB converted
        end
        clear i_trial
    end
    clear ii i_elect tmp_ersp
end
clear i_part



% /////////////////////////////////////////////////////////////////////////
%% Save Data
chanlocs = EEG.chanlocs; %going to want to save electrode locations
% this is large file so it will take some time to save
save([saveLocation 'all_ersp_R_v4.mat'],'all_erspR','chanlocs','-v7.3')

% /////////////////////////////////////////////////////////////////////////
%% OR Load Processed Data If Exists
all_erspR = struct2cell(load([saveLocation 'all_ersp_R_v4.mat'],'all_erspR'));  %gets loaded as a struct
all_erspR = all_erspR{1};
chanlocs = struct2cell(load([saveLocation 'all_ersp_R_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};

% /////////////////////////////////////////////////////////////////////////




% #########################################################################
% /////////////////////////////////////////////////////////////////////////
%% Useful Plots
% /////////////////////////////////////////////////////////////////////////
% #########################################################################



% /////////////////////////////////////////////////////////////////////////
%% # ERS: Power by Model SD #
% /////////////////////////////////////////////////////////////////////////

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Load previously created data (if it has been created)
load([saveLocation 'pwrR_AvG_v4.mat']);
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
        part_ersp = all_erspR{i_part,i_elect}.trials; %get single subject's baseline corrected power
        
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
% #########################################################################
%% Gets a count of trials
err_trl_count(:,1) = cellfun(@numel,x_errdeg_m); %small errors
err_trl_count(:,2) = cellfun(@numel,n_errdeg_m); %large errors
% err_trl_count(:,3) = cell2mat({ALLEEG(1:end).trials}); %total trial count
% ########################################################################
% /////////////////////////////////////////////////////////////////////////
%% Save data if not saved yet
save([saveLocation 'pwrR_AvG_v4.mat'],'errlims','n_errdeg_m','n_pwr','x_errdeg_m','x_pwr');
% /////////////////////////////////////////////////////////////////////////



% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% Plot spectogram across subjects &&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% Raw ERS plots
cmap = redblue(256); %create colormap colors
savename = 'SpecPlotR_';
for ii = 1:length(exp.singletrialselecs)
    
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
    %mean across subjects
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %small errors
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %large errors
    
    CLim = [4.5 6.5]; %set power scale of plot
    
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
    set(get(t,'ylabel'),'String', 'Log Power');
    
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
    set(get(t,'ylabel'),'String', 'Log Power');
    
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
    
    CLim = [-0.1 0.1]; %set power scale of plot
    
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
    set(get(t,'ylabel'),'String', 'Log Power');
    
    savefig([saveFig savename exp.singtrlelec_name{ii}])
    
    clear plot_ers_x plot_ers_n CLim t
end
clear ii i_elect cmap savename


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Grand average difference ERS plot
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


% Raw ERS plots
cmap = jet; %create colormap colors
CLim = [4.5 6.5]; %set power scale of plot
    
% Plot Small Errors
figure('Position', [1 1 1685 405]); colormap(cmap) %open a new figure
subplot(1,2,1)
imagesc(times,freqs,plot_x_avg,CLim);
title(['Accurate: Grand Avg']); set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
% xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Log Power');

% Plot Large Errors
subplot(1,2,2)
imagesc(times,freqs,plot_n_avg,CLim);
title(['Guesses: Grand Avg']); set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
% xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Log Power');

savefig([saveFig 'SpecPlotR_GrandAvg'])

clear cmap savename CLim t


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Plot Accurate-Guesses
CLim = [-0.1 0.1]; %set power scale of plot
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
set(get(t,'ylabel'),'String', 'Log Power');

savefig([saveFig 'SpecPlot_DifR_GrandAvg'])

clear plot_x_avg plot_n_avg CLim t cmap

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------








