% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                               INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% Need the phase opposition code found here: 
%           https://github.com/rufinv/phase-opposition-code
%
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

% And load electrode locations
chanlocs = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};

%% Location to save figures
saveFig = [pwd '\Figures\PhaseOpp\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

%% Load previously processed target-aligned epoch data
% These are large data files and take time to load. Recommend using the
% data files that have been processed which should be on osf, github, or
% can be obtained by request from Sarah Sheldon @ ssheldon@ualberta.ca

all_ersp_trials = load([saveLocation 'all_ersp_v4.mat']); %gets loaded as a struct
all_ersp = all_ersp_trials.all_ersp;
clear all_ersp_trials

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% exp is a function that might get used below
exp2 = exp;
clear exp


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Load phase values separated by trial type using SD (accurate & guess)
load([saveLocation 'phase_trial_v4.mat']); %phase by trial type

% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% OR separate phase by trial type using SD (accurate & guess)
current_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
x_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
n_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
errlims = cell(1,length(exp2.participants));    %pre-allocate
for i_part = 1:length(exp2.participants)
   
    % Get upper and lower limits based on model fit
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp{i_part,i_elect}; %get single subject's ersp
        current_phase{i_part,i_elect}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));

        
        % Get phase values separated by errors in trial
        % small errors
        x_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] );
        % larger errors
        n_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] );
        
        clear part_ersp i_trial
    end
    clear ii
end
clear i_elect i_part

% -------------------------------------------------------------------------
%% Save phase values separated by trial type using SD (accurate & guess)
save([saveLocation 'phase_trial_v4.mat'],'n_phase','x_phase'); %phase by trial type


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load previously processed phase opposition results
load([saveLocation 'POS_pvals_v4_1e5perm.mat']);

% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% OR calculate phase opposition by subject
% -------------------------------------------------------------------------
nperm = 1e5; %number of permutations for non-parametric (POS) test (usually 1000)
for i_part = 1:length(exp2.participants)
%     for ii = 1
    for ii = 1:length(exp2.singletrialselecs)    
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get POS p-vals at each freq x time across trial
        [p_xn.circWW{i_part,i_elect}, p_xn.POS{i_part,i_elect}, p_xn.zPOS{i_part,i_elect}] =...
            PhaseOpposition(x_phase{i_part,i_elect},n_phase{i_part,i_elect},nperm);
        
       clear i_elect
    end
    clear ii
end
clear i_part nperm

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% Combine phase opposition pvals across subjects
% -------------------------------------------------------------------------
% circWW pvals (this is measure used in paper)
% -> better in asymmetric situations where either the relative trial number or 
% the ERP amplitude differed markedly between the two trial groups (VanRullen 2016)
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_xn(i_part,:,:) = p_xn.circWW{i_part,i_elect};
    end
    
    circWWcat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors
    
    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect

% .........................................................................
% POS pvals
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes

    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_xn(i_part,:,:) = p_xn.POS{i_part,i_elect};
    end
    
    POScat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors
    
    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect

% .........................................................................
% zPOS pvals
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_xn(i_part,:,:) = p_xn.zPOS{i_part,i_elect};
    end
    
    zPOScat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors

    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Save phase opposition results
save([saveLocation 'POS_pvals_v4_1e5perm.mat'],'p_xn','zPOScat','POScat',...
    'circWWcat');


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Accurate vs Guess Plots
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% Plots of circWW values at electrodes
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80);
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    figure; colormap(cmap); %open a new figure
    CLim = [0 0.05]; %cuz plotting p-vals

    tmp_plot = squeeze(circWWcat{1,1}(i_elect,:,:)); %extract data
    [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'pdep','yes'); %multi-comp correction
    imagesc(times,freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Accurate v Guesses: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    
    savefig([saveFig 'CircWW_AvG_' exp2.elec_names{ii}])

   clear i_elect i_cond tmp_plot adj_ci_cvrg crit_p h
end
clear ii ncond conds CLim cmap

% -------------------------------------------------------------------------
% Calculate average across electrodes of circWW 
tmp_out = circWWcat{1,1}(2:32,:,:);
circWWcat_avg = combine_pvalues(tmp_out,1); %grand average
clear tmp_out

% Plot grand average
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80);

figure; colormap(cmap); %open a new figure
CLim = [0 0.05]; %cuz plotting p-vals

[h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(circWWcat_avg,0.05,'pdep'); %multi-comp correction
imagesc(times,freqs,tmp_plot,CLim)
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
% xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
set(get(t,'ylabel'),'String', 'p-value');
title(['Accurate v Guesses: Grand Avg'],'FontSize',10.5);

savefig([saveFig 'CircWW_AvG_GrandAvg'])

clear i_elect i_cond tmp_plot adj_ci_cvrg crit_p h CLim cmap


% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
%% Plots of zPOS values at electrodes
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80);
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    figure; colormap(cmap); %open a new figure
    CLim = [0 0.05]; %cuz plotting p-vals

    tmp_plot = squeeze(zPOScat{1,1}(i_elect,:,:)); %extract data
    [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'dep','yes'); %multi-comp correction
    imagesc(times,freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Accurate v Guesses: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    
    savefig([saveFig 'zPOS_AvG_' exp2.elec_names{ii}])

    clear i_elect i_cond tmp_plot adj_ci_cvrg crit_p h
end
clear ii ncond conds CLim cmap

% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
%% Plots of POS values at electrodes
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80);
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    figure; colormap(cmap); %open a new figure
    CLim = [0 0.05]; %cuz plotting p-vals

    tmp_plot = squeeze(POScat{1,1}(i_elect,:,:)); %extract data
    [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'pdep','yes'); %multi-comp correction
    imagesc(times,freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([2 40]); yticks(5:5:40)
%     xlim([-200 800]); xticks(-200:100:800) %match ERPs
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Accurate v Guesses: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    
    savefig([saveFig 'POS_AvG_' exp2.elec_names{ii}])

   clear i_elect i_cond tmp_plot adj_ci_cvrg crit_p h
end
clear ii ncond conds CLim cmap
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
% freqband = [30 40]; %gamma
% freqband = [23 29]; %beta2
% freqband = [15 22]; %beta1
% freqband = [8 14]; %alpha
freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));


% Choose phase measure
cout_pmap = circWWcat{1,1};

% Get mean power at frequency band for each electrode
pwr_top = NaN([length(times),length(elect_erp),1]); %pre-allocate
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    tmp_plot = squeeze(cout_pmap(i_elect,:,:)); %extract data
    [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'pdep','yes'); %multi-comp correction
    
    pwr_top(:,i_elect,1) = squeeze(mean(tmp_plot(freqlim,:),1));
    clear tmp_plot h crit_p adj_ci_cvrg
end
clear ii i_elect

figname = 'circWW_TopPlot_'; %name of figure
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],60); %create colormap colors
CLim = [0 0.5]; %color axis limits

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
    set(get(t,'ylabel'),'String', 'p-value');

    savefig([saveFig figname num2str(freqband(1)) '-' num2str(freqband(2)) '_' num2str(itWin(1)) 'to' num2str(itWin(2))])
    
    clear itWin time_window temp
end
clear tw_i t cmap

clear freqlim freqband CLim pwr_top cout_pmap


clear tWin elect_erp figname 




% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% OPTIONAL: Load previously processed catch trial data
% Catch trials were not included in the phase analysis of the paper, but
% can be explored using this code.
% Can use LoadProcData_OrientTask.m function to directly load the data if
% there is no saved file yet
% 
% assumes that the catch trials have been processed with the same
% parameters as the all_ersp data (check their settings in exp to make sure) 
catch_trials_v1_wav = load('all_ersp_catch_trials_v2.mat');
all_ersp_byC = catch_trials_v1_wav.all_ersp;

% Check to make sure their frequency and time scales are the same
% output should be 0 or they are not the same!
sum(catch_trials_v1_wav.freqs ~= freqs)
sum(catch_trials_v1_wav.times ~= times)
clear catch_trials_v1_wav

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Load phase values separated by trial type using SD (accurate & guess)
load([saveLocation 'phase_trial_v4.mat']); %phase by trial type

% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% OR separate phase by trial type using SD (accurate & guess)
% -------------------------------------------------------------------------
current_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
x_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
n_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
errlims = cell(1,length(exp2.participants));    %pre-allocate
for i_part = 1:length(exp2.participants)
   
    % Get upper and lower limits based on model fit
    errlims{i_part}(1) = -(model_out{1,i_part}(2)); %negative value
    errlims{i_part}(2) = model_out{1,i_part}(2);
    
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp{i_part,i_elect}; %get single subject's ersp
        current_phase{i_part,i_elect}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));

        
        % Get phase values separated by errors in trial
        % small errors
        x_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find((resp_errdeg{i_part}<(errlims{i_part}(2)*0.75) & resp_errdeg{i_part}>(errlims{i_part}(1)*0.75)))] );
        % larger errors
        n_phase{i_part,i_elect}(:,:,:) = current_phase{i_part,i_elect}(:,:,...
            [find(resp_errdeg{i_part}>=(errlims{i_part}(2)*1.5)) find(resp_errdeg{i_part}<=(errlims{i_part}(1)*1.5))] );
        
        clear part_ersp i_trial
    end
    clear ii
end
clear i_elect i_part

% .............
% Phase for catch trials
c_phase = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp2.participants)
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp_byC{i_part,i_elect}; %get single subject's ersp

        % Get phase values from catch trials
        c_phase{i_part,i_elect}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));
        
        clear part_ersp
    end
    clear ii
end
clear i_elect i_part

% -------------------------------------------------------------------------
%% Save phase values separated by trial type using SD (accurate & guess)
save([saveLocation 'phase_trial_v4.mat'],'n_phase','x_phase','c_phase'); %phase by trial type

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Calculate p-value by subject
% -------------------------------------------------------------------------
nperm = 1e5; %number of permutations for non-parametric (POS) test (usually 1000)
for i_part = 1:length(exp2.participants)
%     for ii = 1
    for ii = 1:length(exp2.singletrialselecs)    
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get POS p-vals at each freq x time across trial
        [p_xn.circWW{i_part,i_elect}, p_xn.POS{i_part,i_elect}, p_xn.zPOS{i_part,i_elect}] =...
            PhaseOpposition(x_phase{i_part,i_elect},n_phase{i_part,i_elect},nperm);
        [p_cn.circWW{i_part,i_elect}, p_cn.POS{i_part,i_elect}, p_cn.zPOS{i_part,i_elect}] =...
            PhaseOpposition(c_phase{i_part,i_elect},n_phase{i_part,i_elect},nperm);
        [p_cx.circWW{i_part,i_elect}, p_cx.POS{i_part,i_elect}, p_cx.zPOS{i_part,i_elect}] =...
            PhaseOpposition(c_phase{i_part,i_elect},x_phase{i_part,i_elect},nperm);
        
       clear i_elect
    end
    clear ii
end
clear i_part nperm

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% Combine phase opposition pvals across subjects
% -------------------------------------------------------------------------
% circWW pvals
% -> better in asymmetric situations where either the relative trial number or 
% the ERP amplitude differed markedly between the two trial groups (VanRullen 2016)
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    c_pval_cn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_cx = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_cn(i_part,:,:) = p_cn.circWW{i_part,i_elect};
        c_pval_cx(i_part,:,:) = p_cx.circWW{i_part,i_elect};
        c_pval_xn(i_part,:,:) = p_xn.circWW{i_part,i_elect};
    end
    
    circWWcat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors
    circWWcat{2}(i_elect,:,:) = combine_pvalues(c_pval_cx,1); %catch vs small errors
    circWWcat{3}(i_elect,:,:) = combine_pvalues(c_pval_cn,1); %catch vs large errors
    
    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect

% .........................................................................
% POS pvals
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    c_pval_cn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_cx = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_cn(i_part,:,:) = p_cn.POS{i_part,i_elect};
        c_pval_cx(i_part,:,:) = p_cx.POS{i_part,i_elect};
        c_pval_xn(i_part,:,:) = p_xn.POS{i_part,i_elect};
    end
    
    POScat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors
    POScat{2}(i_elect,:,:) = combine_pvalues(c_pval_cx,1); %catch vs small errors
    POScat{3}(i_elect,:,:) = combine_pvalues(c_pval_cn,1); %catch vs large errors
    
    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect

% .........................................................................
% zPOS pvals
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    c_pval_cn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_cx = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    c_pval_xn = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    for i_part = 1:length(exp2.participants)
        c_pval_cn(i_part,:,:) = p_cn.zPOS{i_part,i_elect};
        c_pval_cx(i_part,:,:) = p_cx.zPOS{i_part,i_elect};
        c_pval_xn(i_part,:,:) = p_xn.zPOS{i_part,i_elect};
    end
    
    zPOScat{1}(i_elect,:,:) = combine_pvalues(c_pval_xn,1); %small vs large errors
    zPOScat{2}(i_elect,:,:) = combine_pvalues(c_pval_cx,1); %catch vs small errors
    zPOScat{3}(i_elect,:,:) = combine_pvalues(c_pval_cn,1); %catch vs large errors

    clear c_pval_cn c_pval_cx c_pval_xn i_part   
end
clear ii i_elect


% -------------------------------------------------------------------------
%% Save phase opposition results
save([saveLocation 'POS_pvals_v4_1e5perm.mat'],'p_cn','p_xn','p_cx','POScat',...
    'circWWcat','zPOScat');


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------


clear saveFig








