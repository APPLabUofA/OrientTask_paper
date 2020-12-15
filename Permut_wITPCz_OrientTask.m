% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                               INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% This is calculates the weighted ITPC and standardized the data with the
% null hypothesis distribution to get statistical meaningful results.
% Further details of the method can be found in Cohen's book "Analyzing 
% neural time series data: theory and practice"
% 
% Need the combine_pvalues function found at link below or can use the
% Stouffer's equation provided in the comments 
%           https://github.com/rufinv/phase-opposition-code
% 
% This method takes a lot of time to run depending on your computer and the
% number of permutations you choose so be aware that this can tie up your 
% computer for days!
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
saveFig = [pwd '\Figures\wITPCz\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% exp is a function that gets used below
exp2 = exp;
clear exp

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Load previously processed data if no phase .mat file has been created
% These are large data files and take time to load. Recommend using the
% data files that have been processed which should be on osf, github, or
% can be obtained by request from Sarah Sheldon at ssheldon@ualberta.ca

all_ersp_trials = load([saveLocation 'all_ersp_v4.mat']); %gets loaded as a struct
all_ersp = all_ersp_trials.all_ersp;
clear all_ersp_trials

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Load phase values if it exists
load([saveLocation 'phase_angle.mat']); 

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% OR extract phase data from 

phase_angle = cell(length(exp2.singletrialselecs),length(exp2.participants)); %pre-allocate
for i_part = 1:length(exp2.participants)
   
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp(:,:,:) = all_ersp{i_part,i_elect}; %get single subject's ersp
        phase_angle{i_elect,i_part}(:,:,:) = squeeze(angle(part_ersp(:,:,:)));

        clear part_ersp i_trial
    end
    clear ii
end
clear i_elect i_part

% -------------------------------------------------------------------------
%% Save phase values
save([saveLocation 'phase_angle.mat'],'phase_angle'); 


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load previously analyzed wITPC_z data
load([saveLocation 'witpcz_out_1e4perm.mat'])



% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Calculate wITPC_z maps per subject
% wITPC is response error modulating the length of phase angles

% specify values  
n_permutes = 1e4; %number of permutations used to estimate null distribution (takes really long time with increasing number)
% voxel_pval = 0.05; %uncorrect pixel-level threshold
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction

% pre-allocate variables
witpc = NaN(length(exp2.participants),length(exp2.singletrialselecs),length(freqs),length(times)); %pre-allocate
zmap_out = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_witpc_thresh = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
witpc_threshold_out = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
threshold_out = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_z = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
% Perform analysis on each participant
for i_part = 1:length(exp2.participants)
   
    tmp_err = abs(resp_errdeg{i_part}); %abs value of response errors
    n_resp = length(tmp_err); %number of trials
    
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        %create a freq x time matrix of the response error values
        tmp_resperr = permute(repmat(tmp_err,[length(times) 1 length(freqs)]),[3 1 2]);
        
        %compute wITPC
        witpc(i_part,ii,:,:) = squeeze(abs(nanmean(tmp_resperr.*squeeze(exp(1i*phase_angle{i_elect,i_part})),3)));
        
        clear tmp_resperr
        
        
        %% Create H0 distribution
        perm_witpc  = NaN(n_permutes,length(freqs),length(times)); % initialize null hypothesis matrices
        max_pixel_vals = NaN(n_permutes,1); % initialize null hypothesis matrices
        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:n_permutes
            fake_err_map = tmp_err(randperm(n_resp)); %randomize order of response errors

            %create a freq x time matrix of the response error values
            tmp_resperr = permute(repmat(fake_err_map,[length(times) 1 length(freqs)]),[3 1 2]);

            % save all permuted values
            perm_witpc(permi,:,:) = squeeze(abs(nanmean(tmp_resperr.*squeeze(exp(1i*phase_angle{i_elect,i_part})),3)));

            max_pixel_vals(permi) = max(max(perm_witpc(permi,:,:)));
            
            clear fake_err_map tmp_resperr
        end
        clear permi
        
        %% Real z-map
        % now compute Z-map (standardized units from each value away from
        % the distribution of null-hypothesis values)
        permmean = squeeze(mean(perm_witpc,1));
        permstd  = squeeze(std(perm_witpc,[],1));
        witpc_z = (squeeze(squeeze(witpc(i_part,ii,:,:)))-permmean)./permstd;
        zmap_out{i_part,ii} = witpc_z; %save values of plot
        
        % Z-value to p %(p_z method)
        p_witpc_z{i_part,ii} = 2*(1-normcdf(abs(witpc_z)));  %times 2 for 2-tailed test
        
        % just wITPC
        upper_threshold = prctile(max_pixel_vals(:),100-mcc_voxel_pval*100);
        witpc_threshold_out{i_part,ii} = upper_threshold;
        zmapthresh = witpc_z;
        zmapthresh(squeeze(squeeze(witpc(i_part,ii,:,:)))<upper_threshold)=0;
        zmap_witpc_thresh{i_part,ii} = zmapthresh; %save values of plot
        
        clear max_pixel_vals upper_threshold zmapthresh
        
        %% Create pixel-level corrected threshold
        max_pixel_vals = NaN(n_permutes,2); % initialize null hypothesis matrices
        for permi = 1:n_permutes
            % null z-map
            fakez = (squeeze(perm_witpc(permi,:,:))-permmean) ./ permstd;
            
            % save maximum pixel values
            max_pixel_vals(permi,:) = [ min(fakez(:)) max(fakez(:)) ];

            clear fakez
        end
        clear permi permmean permstd
        
        % apply pixel-level corrected threshold
        lower_threshold = prctile(max_pixel_vals(:,1),    mcc_voxel_pval*100/2);
        upper_threshold = prctile(max_pixel_vals(:,2),100-mcc_voxel_pval*100/2);
        threshold_out{i_part,ii} = [lower_threshold upper_threshold];
        zmapthresh = witpc_z;
        zmapthresh(witpc_z>lower_threshold & witpc_z<upper_threshold)=0;
       
        clear perm_witpc max_pixel_vals lower_threshold upper_threshold zmapthresh...
            i_elect witpc_z
    end
    clear ii n_resp tmp_err
end
clear i_part voxel_pval mcc_cluster_pval mcc_voxel_pval n_permutes


% -------------------------------------------------------------------------
%% Save analyzed wITPC_z data
% some variables that are saved were not used in the paper
save([saveLocation 'witpcz_out_1e4perm.mat'],'witpc','zmap_out','zmap_witpc_thresh',...
    'witpc_threshold_out','threshold_out','p_witpc_z'); 



% /////////////////////////////////////////////////////////////////////////
%% Combine pvals across subjects
% /////////////////////////////////////////////////////////////////////////
%Stouffer
% result = squeeze(1-normcdf(sum(norminv(1-pmatrix),dim)./sqrt(size(pmatrix,dim))));
% .........................................................................

tmp_pcorr = p_witpc_z; %pvals to combine

pcorr_cat = NaN([length(exp2.singletrialselecs),length(freqs),length(times)]); %pre-allocate
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    combo_p_z = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    combo_z = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate

    for i_part = 1:length(exp2.participants)
        combo_p_z(i_part,:,:) = tmp_pcorr{i_part,ii};
        combo_z(i_part,:,:) = zmap_out{i_part,ii};
    end
    
    pcorr_cat(i_elect,:,:) = combine_pvalues(combo_p_z,1,1);

    clear combo_p_z i_part combo_z   
end
clear ii i_elect tmp_pcorr

% -------------------------------------------------------------------------
% Calculate average across electrodes
tmp_out = pcorr_cat(2:32,:,:);
pcorr_cat_avg = squeeze(mean(tmp_out,1)); %grand average
clear tmp_out


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Plots of combined p-values at electrodes
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80); %white to green
tite = 'p_wITPCz: ';
savname = 'p_wITPCz_';
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    figure; colormap(cmap); %open a new figure
    CLim = [0 0.05]; %cuz plotting p-vals

    tmp_plot = squeeze(pcorr_cat(i_elect,:,:)); %extract data
    [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'dep','yes'); %multi-comp correction
    tmp_plot(tmp_plot>0.05)=1; %set non-significant values to 1
    imagesc(times,freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
%     ylim([3 40]); yticks(5:5:40)
%     xlim([-700 800]); xticks(-600:200:800)
    ylim([2 40]); yticks(5:5:40)
    xlim([-200 800]); xticks(-200:100:800) %match ERPs
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title([tite exp2.singtrlelec_name{ii}],'FontSize',10.5);
    
    savefig([saveFig savname exp2.elec_names{ii}])

   clear i_elect i_cond tmp_plot adj_ci_cvrg crit_p h
end
clear ii ncond conds CLim cmap tite t savname


% /////////////////////////////////////////////////////////////////////////
% Plot Grand Average
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80); %used in paper
figure; colormap(cmap); %open a new figure
CLim = [0 0.05]; %cuz plotting p-vals
tite = 'p_wITPCz: ';
savname = 'p_wITPCz_';

tmp_plot = pcorr_cat_avg; %extract data
[h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'dep','yes'); %multi-comp correction
tmp_plot(tmp_plot>0.05)=1; %set non-significant values to 1
imagesc(times,freqs,tmp_plot,CLim)
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
% ylim([3 40]); yticks(5:5:40)
% xlim([-700 800]); xticks(-600:200:800)
ylim([2 40]); yticks(5:5:40)
xlim([-200 800]); xticks(-200:100:800) %match ERPs
ylabel('Freqency (Hz)'); xlabel('Time (ms)');
t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
set(get(t,'ylabel'),'String', 'p-value');
title([tite 'Grand Avg'],'FontSize',10.5);

savefig([saveFig savname 'GrandAvg'])

clear tmp_plot adj_ci_cvrg crit_p h CLim cmap tite t savname




% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     ''''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Not used in paper, but good for visualization
% .........................................................................

tmp_pcorr = p_witpc_z; %pvals to combine

pcorr_cat = NaN([length(exp2.singletrialselecs),length(freqs),length(times)]); %pre-allocate
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    combo_p_z = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate
    combo_z = NaN([length(exp2.participants),length(freqs),length(times)]); %pre-allocate

    for i_part = 1:length(exp2.participants)
        combo_p_z(i_part,:,:) = tmp_pcorr{i_part,ii};
        combo_z(i_part,:,:) = zmap_out{i_part,ii};
    end
    
    pcorr_cat(i_elect,:,:) = combine_pvalues(combo_p_z,1,1);
    
    clear combo_p_z i_part combo_z   
end
clear ii i_elect tmp_pcorr

% -------------------------------------------------------------------------

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
freqband = [23 29]; %beta2
% freqband = [15 22]; %beta1
% freqband = [10 14]; %high alpha
% freqband = [8 14]; %alpha
% freqband = [8 11]; %low alpha
% freqband = [4 7]; %theta
% freqband = [2 3]; %delta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));


cout_pmap = pcorr_cat;
%note that this fdr multiple corrections is more conservative than the plots above
% [h, crit_p, adj_ci_cvrg, cout_pmap] = fdr_bh(cout_pmap,0.05,'pdep','yes'); %multi-comp correction
% clear h crit_p adj_ci_cvrg

% Get mean power at frequency band for each electrode
pwr_top = NaN([length(times),length(elect_erp),1]); %pre-allocate
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    pwr_top(:,i_elect,1) = squeeze(mean(cout_pmap(i_elect,freqlim,:),2));
end
clear ii i_elect

% Settings for plot
figname = 'p_wITPCz_Topo_'; %name of figure
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],60); %create colormap colors
CLim = [0 0.05]; %color axis limits

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
 

% -------------------------------------------------------------------------

% clears variables that end/begin with...
clear -regexp \<p_ \<zmap_ _thresh\>

clear cout_thresh cout_thresh pcorr_cat cout_zmap pcorr_cat_avg saveFig...
    witpc witpc_threshold_out threshold_out el_erp_names

% -------------------------------------------------------------------------













