% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% Need to have run Mdn_PwrR_Split.m to use the code below.
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
load([saveLocation 'mdn_pwrR_split_v4.mat'])

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
%%    Compute mean g & sd in frequency windows then plot
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% Going to need to uncomment one frequency band at a time and then run
% the code below it (I know it is a bit annoying, but it was to make sure
% visual checks happened throughout the analysis)

clear g_Hpwr_win g_Lpwr_win sd_Hpwr_win sd_Lpwr_win freqlim freqband...
    g_out_bypwr sd_out_bypwr perm_sd_pwr_win perm_g_pwr_win

%finds the frequencies you want (gamma (30–90 Hz))
% freqband = [30 40]; %gamma
% freqband = [15 29]; %beta
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

g_Hpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
g_Lpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
sd_Hpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
sd_Lpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % (participant x electrode x freq x time)
        g_Hpwr_win(i_part,i_elect,:) = squeeze(mean(g_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        g_Lpwr_win(i_part,i_elect,:) = squeeze(mean(g_out_Lpwr(i_part,i_elect,freqlim,:),3));
        sd_Hpwr_win(i_part,i_elect,:) = squeeze(mean(sd_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        sd_Lpwr_win(i_part,i_elect,:) = squeeze(mean(sd_out_Lpwr(i_part,i_elect,freqlim,:),3));
    end
    clear ii i_elect
end
clear i_part

% -------------------------------------------------------------------------
% Average g across subjects
g_out_bypwr(:,:,1) = squeeze(mean(g_Hpwr_win(:,:,:),1)); 
g_out_bypwr(:,:,2) = squeeze(mean(g_Lpwr_win(:,:,:),1)); 
g_out_bypwr(:,:,3) = squeeze(mean((g_Hpwr_win(:,:,:)-g_Lpwr_win(:,:,:)),1)); %difference
% .........................................
% Average sd across subjects
sd_out_bypwr(:,:,1) = squeeze(mean(sd_Hpwr_win(:,:,:),1)); 
sd_out_bypwr(:,:,2) = squeeze(mean(sd_Lpwr_win(:,:,:),1)); 
sd_out_bypwr(:,:,3) = squeeze(mean((sd_Hpwr_win(:,:,:)-sd_Lpwr_win(:,:,:)),1)); %difference
% ------------------------------------------------------------------------- 

% Plot g with error bars
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    % get axes limits
    ymin = -0.1; ymax = 0.3;
    xmin = -700; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(times,g_out_bypwr(i_elect,:,1),squeeze(std(g_Hpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            times,g_out_bypwr(i_elect,:,2),squeeze(std(g_Lpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m',...
            times,g_out_bypwr(i_elect,:,3),squeeze(std((g_Hpwr_win(:,i_elect,:)-g_Lpwr_win(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'k');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    xlim([xmin xmax]); xticks(-600:200:800)
    
    title([exp.elec_names{ii} ': g ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz']); 
    xlabel('Time (ms)'); ylabel('Guess Rate (g)'); 
        
    legend({'High Power','Low Power','Difference'},'Location','best');
    
    savefig([saveFigLoc 'g_' exp.elec_names{ii} '_' num2str(freqband(1)) '-' num2str(freqband(2))])

    hold off
end
clear ii xmax xmin ymin ymax i_elect 
clear g_out_bypwr

% /////////////////////////////////////////////////////////////////////////
% Plot grand average g with error bars

% .........................................
% Grand average g across subjects & electrodes
g_out_bypwr(:,1) = squeeze(mean(mean(g_Hpwr_win(:,2:32,:),2),1)); 
g_out_bypwr(:,2) = squeeze(mean(mean(g_Lpwr_win(:,2:32,:),2),1)); 
g_out_bypwr(:,3) = squeeze(mean(mean((g_Hpwr_win(:,2:32,:)-g_Lpwr_win(:,2:32,:)),2),1)); %difference
% .........................................

% get axes limits
ymin = -0.1; ymax = 0.3;
xmin = -700; xmax = 800;

figure('Color',[1 1 1]); 
boundedline(times,g_out_bypwr(:,1),squeeze(std(squeeze(mean(g_Hpwr_win(:,2:32,:),2)),[],1))./sqrt(length(exp.participants)),'c',...
        times,g_out_bypwr(:,2),squeeze(std(squeeze(mean(g_Lpwr_win(:,2:32,:),2)),[],1))./sqrt(length(exp.participants)),'m',...
        times,g_out_bypwr(:,3),squeeze(std(squeeze(mean((g_Hpwr_win(:,2:32,:)-g_Lpwr_win(:,2:32,:)),2)),[],1))./sqrt(length(exp.participants)),'k');
hold on
line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
xlim([xmin xmax]); xticks(-600:200:800)

title(['Grand Avg: g ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz']); 
xlabel('Time (ms)'); ylabel('Guess Rate (g)'); 

legend({'High Power','Low Power','Difference'},'Location','best');

savefig([saveFigLoc 'g_GrandAvg_' num2str(freqband(1)) '-' num2str(freqband(2))])

hold off

clear ii xmax xmin ymin ymax i_elect 
clear g_out_bypwr


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////

% Plot sd with error bars
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    % get axes limits
    ymin = -2; ymax = 14;
    xmin = -700; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(times,sd_out_bypwr(i_elect,:,1),squeeze(std(sd_Hpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            times,sd_out_bypwr(i_elect,:,2),squeeze(std(sd_Lpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m',...
            times,sd_out_bypwr(i_elect,:,3),squeeze(std((sd_Hpwr_win(:,i_elect,:)-sd_Lpwr_win(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'k');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    xlim([xmin xmax]); xticks(-600:200:800)
    
    title([exp.elec_names{ii} ': SD ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz']); 
    xlabel('Time (ms)'); ylabel('Variability (sd)'); 
        
    legend({'High Power','Low Power','Difference'},'Location','best');
    
    savefig([saveFigLoc 'sd_' exp.elec_names{ii} '_' num2str(freqband(1)) '-' num2str(freqband(2))])

    hold off
end
clear ii xmax xmin ymin ymax i_elect



% ------------------------------------------------------------------------- 
clear g_Hpwr_win g_Lpwr_win sd_Hpwr_win sd_Lpwr_win freqlim freqband...
    g_out_bypwr sd_out_bypwr
% -------------------------------------------------------------------------




% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% :::::::::::::::::::::  Permutation Test  :::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 
% **Need code from: https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox

% ------------------------------------------------------------------------- 
%% Guess rate

clear g_Hpwr_win g_Lpwr_win perm_g_pwr_win g_out_bypwr freqlim freqband

%finds the frequencies you want (gamma (30–90 Hz))
freqband = [30 40]; %gamma
% freqband = [15 29]; %beta
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

g_Hpwr_win = NaN(length(exp.singletrialselecs),length(times),length(exp.participants)); %pre-allocate
g_Lpwr_win = NaN(length(exp.singletrialselecs),length(times),length(exp.participants)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
                                       % (participant x electrode x freq x time)
        g_Hpwr_win(i_elect,:,i_part) = squeeze(mean(g_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        g_Lpwr_win(i_elect,:,i_part) = squeeze(mean(g_out_Lpwr(i_part,i_elect,freqlim,:),3));
    end
    clear ii i_elect
end
clear i_part

% ------------------------------------------------------------------------- 
% Test only relevant times
timewin = [-700 800];
timelim = find(times>=timewin(1) & times<=timewin(2));
% ------------------------------------------------------------------------- 

% One sample/repeated-measures permutation test
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(...
    g_Hpwr_win(2:32,timelim,:)-g_Lpwr_win(2:32,timelim,:)),1e5,.01);

% Save results where figures are for later plotting
mat_name = ['g_mdnPwr_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
save([saveFigLoc mat_name],'pval','t_orig','tmx_ptile')

clear pval t_orig tmx_ptile seed_state est_alpha mat_name


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% One sample/repeated-measures permutation test - Grand Average
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(mean((...
    g_Hpwr_win(2:32,timelim,:)-g_Lpwr_win(2:32,timelim,:)),1),1e5,.01);

% Save results where figures are for later plotting
mat_name = ['g_mdnPwr_GrandAvg_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
save([saveFigLoc mat_name],'pval','t_orig','tmx_ptile')

clear pval t_orig tmx_ptile seed_state est_alpha mat_name


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% clears variables that end/begin with...
clear -regexp \<permtest_ \<ttest_
clear g_Hpwr_win g_Lpwr_win perm_g_pwr_win g_out_bypwr freqlim freqband...
    timewin timelim

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ------------------------------------------------------------------------- 
%% SD parameter

clear sd_Hpwr_win sd_Lpwr_win freqlim freqband sd_out_bypwr perm_sd_pwr_win

%finds the frequencies you want (gamma (30–90 Hz))
freqband = [30 40]; %gamma
% freqband = [15 29]; %beta
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

sd_Hpwr_win = NaN(length(exp.singletrialselecs),length(times),length(exp.participants)); %pre-allocate
sd_Lpwr_win = NaN(length(exp.singletrialselecs),length(times),length(exp.participants)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
                                       % (participant x electrode x freq x time)
        sd_Hpwr_win(i_elect,:,i_part) = squeeze(mean(sd_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        sd_Lpwr_win(i_elect,:,i_part) = squeeze(mean(sd_out_Lpwr(i_part,i_elect,freqlim,:),3));
    end
    clear ii i_elect
end
clear i_part

% ------------------------------------------------------------------------- 
% Test only relevant times
timewin = [-700 800];
timelim = find(times>=timewin(1) & times<=timewin(2));
% ------------------------------------------------------------------------- 

% One sample/repeated-measures permutation test
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(...
    sd_Hpwr_win(2:32,timelim,:)-sd_Lpwr_win(2:32,timelim,:)),1e5,0.01);

% Save results where figures are for later plotting
mat_name = ['sd_mdnPwr_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
save([saveFigLoc mat_name],'pval','t_orig','tmx_ptile')

clear pval t_orig tmx_ptile seed_state est_alpha mat_name


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% One sample/repeated-measures permutation test - Grand Average
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(mean((...
    sd_Hpwr_win(2:32,timelim,:)-sd_Lpwr_win(2:32,timelim,:)),1),1e5,.01);

% Save results where figures are for later plotting
mat_name = ['sd_mdnPwr_GrandAvg_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
save([saveFigLoc mat_name],'pval','t_orig','tmx_ptile')

clear pval t_orig tmx_ptile seed_state est_alpha mat_name


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% clears variables that end/begin with...
clear -regexp \<permtest_ \<ttest_
clear sd_Hpwr_win sd_Lpwr_win freqlim freqband sd_out_bypwr...
    perm_sd_pwr_win timewin timelim

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% ::::::::::::::::::  Add Significant 2 Plots  :::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 

%% Guess Rate
%finds the frequencies you want (gamma (30–90 Hz))
% freqband = [30 40]; %gamma
% freqband = [15 29]; %beta
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

g_Hpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
g_Lpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % (participant x electrode x freq x time)
        g_Hpwr_win(i_part,i_elect,:) = squeeze(mean(g_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        g_Lpwr_win(i_part,i_elect,:) = squeeze(mean(g_out_Lpwr(i_part,i_elect,freqlim,:),3));
    end
    clear ii i_elect
end
clear i_part

% -------------------------------------------------------------------------
% Average g across subjects by errors
g_out_bypwr(:,:,1) = squeeze(mean(g_Hpwr_win(:,:,:),1)); 
g_out_bypwr(:,:,2) = squeeze(mean(g_Lpwr_win(:,:,:),1)); 
g_out_bypwr(:,:,3) = squeeze(mean((g_Hpwr_win(:,:,:)-g_Lpwr_win(:,:,:)),1)); %difference
% ------------------------------------------------------------------------- 

% Load previous results from permutation test
mat_name = ['g_mdnPwr_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
load([saveFigLoc mat_name],'pval')

% Find relevant times in full epoch
timewin = [-700 800];
timelim = find(times>=timewin(1) & times<=timewin(2));

% ------------------------------------------------------------------------- 
% Plot g with error bars
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    % get axes limits
    ymin = -0.1; ymax = 0.3;
    xmin = -700; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(times,g_out_bypwr(i_elect,:,1),squeeze(std(g_Hpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            times,g_out_bypwr(i_elect,:,2),squeeze(std(g_Lpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m',...
            times,g_out_bypwr(i_elect,:,3),squeeze(std((g_Hpwr_win(:,i_elect,:)-g_Lpwr_win(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'k');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    xlim([xmin xmax]); xticks(-600:200:800)
    
    title([exp.elec_names{ii} ': g ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz']); 
    xlabel('Time (ms)'); ylabel('Guess Rate (g)'); 
        
    legend({'High Power','Low Power','Difference'},'Location','best');
    
    % Select time points with significant difference
    sigloc = timelim(find(pval(ii,:)<0.05)); %find time points of sig
    gcf; plot(times(sigloc),g_out_bypwr(i_elect,sigloc,3),'.r') %mark points on plot
    clear sigloc
    
    savefig([saveFigLoc 'g_sig_' exp.elec_names{ii} '_' num2str(freqband(1)) '-' num2str(freqband(2))])

    hold off
end
clear ii xmax xmin ymin ymax i_elect 


clear pval t_orig tmx_ptile seed_state est_alpha mat_name timewin timelim...
    g_out_bypwr g_Hpwr_win g_Lpwr_win freqband freqlim



% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
%% SD

%finds the frequencies you want (gamma (30–90 Hz))
% freqband = [30 40]; %gamma
% freqband = [15 29]; %beta
% freqband = [8 14]; %alpha
% freqband = [4 7]; %theta
freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));

sd_Hpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
sd_Lpwr_win = NaN(length(exp.participants),length(exp.singletrialselecs),length(times)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % (participant x electrode x freq x time)
        sd_Hpwr_win(i_part,i_elect,:) = squeeze(mean(sd_out_Hpwr(i_part,i_elect,freqlim,:),3)); 
        sd_Lpwr_win(i_part,i_elect,:) = squeeze(mean(sd_out_Lpwr(i_part,i_elect,freqlim,:),3));
    end
    clear ii i_elect
end
clear i_part

% -------------------------------------------------------------------------
% Average g across subjects by errors
sd_out_bypwr(:,:,1) = squeeze(mean(sd_Hpwr_win(:,:,:),1)); 
sd_out_bypwr(:,:,2) = squeeze(mean(sd_Lpwr_win(:,:,:),1)); 
sd_out_bypwr(:,:,3) = squeeze(mean((sd_Hpwr_win(:,:,:)-sd_Lpwr_win(:,:,:)),1)); %difference
% ------------------------------------------------------------------------- 

% Load previous results from permutation test
mat_name = ['sd_mdnPwr_' num2str(freqband(1)) '-' num2str(freqband(2)) 'Hz_permtest.mat'];
load([saveFigLoc mat_name],'pval')

% Find relevant times in full epoch
timewin = [-700 800];
timelim = find(times>=timewin(1) & times<=timewin(2));

% ------------------------------------------------------------------------- 
% Plot sd with error bars
for ii = 1:length(exp.singletrialselecs)
    i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

    % get axes limits
    ymin = -4; ymax = 14;
    xmin = -700; xmax = 800;
    
    figure('Color',[1 1 1]); 
    boundedline(times,sd_out_bypwr(i_elect,:,1),squeeze(std(sd_Hpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            times,sd_out_bypwr(i_elect,:,2),squeeze(std(sd_Lpwr_win(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m',...
            times,sd_out_bypwr(i_elect,:,3),squeeze(std((sd_Hpwr_win(:,i_elect,:)-sd_Lpwr_win(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'k');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    xlim([xmin xmax]); xticks(-600:200:800)
    
    title([exp.elec_names{ii} ': SD ' num2str(freqband(1)) '-' num2str(freqband(2)) ' Hz']); 
    xlabel('Time (ms)'); ylabel('Standard Deviation (SD)'); 
        
    legend({'High Power','Low Power','Difference'},'Location','best');
    
    % Select time points with significant difference
    sigloc = timelim(find(pval(ii,:)<0.05)); %find time points of sig
    gcf; plot(times(sigloc),sd_out_bypwr(i_elect,sigloc,3),'.r') %mark points on plot
    clear sigloc
    
    savefig([saveFigLoc 'sd_sig_' exp.elec_names{ii} '_' num2str(freqband(1)) '-' num2str(freqband(2))])

    hold off
end
clear ii xmax xmin ymin ymax i_elect 


clear pval t_orig tmx_ptile seed_state est_alpha mat_name timewin timelim...
    sd_out_bypwr sd_Hpwr_win sd_Lpwr_win freqband freqlim



