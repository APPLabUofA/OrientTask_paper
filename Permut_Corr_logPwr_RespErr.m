% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% Does multiple runs of the nonparametric permutation version of Spearman's
% correlation between individual raw power and response errors across the
% time-frequency space. Uses pixel-based multiple corrections procedure.
% 
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Location to save figures
saveFig = [pwd '\Figures\Corr_logPwr_RespErr\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

%% Load previously processed target-aligned epoch data
%data created in Power_RawLog_ERS_Plots.m
all_erspR = struct2cell(load('all_ersp_R_v4.mat'));  %gets loaded as a struct
ERSP = all_erspR{1};

clear all_erspR all_ersp_Z

chanlocs = struct2cell(load([saveLocation 'all_ersp_R_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};


%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Permutation Testing
%set parameters inside for loop

% voxel_pval = 0.05; %uncorrect pixel-level threshold
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction

n_permutes = 1e5; %number of permutations used to estimate null distribution

%finds the times you want from the times variable
timewin = [-700 800];
timephi = find(times>=timewin(1) & times<=timewin(2));

num_frex = length(freqs); %number of frequencies
nTimepoints = length(timephi); %number of timepoints


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

zmap_out = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
threshold_out = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
realcorrs_out = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
realcorrs_plot_out = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
tvalsRho = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
p_zcorr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
p_zcorr_thresh = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
    
    % Get rank-transform response errors
    % errdeg_rank = tiedrank(resp_errdeg{i_part});
    errdeg_rank = tiedrank(abs(resp_errdeg{i_part})); %assume errors to the left = errors to the right
    
    n_resp = length(errdeg_rank); %total number of trials for subject
    
    for ii = 1:length(exp.singletrialselecs)

        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        % Get subject's power
              %all_ersp(participant x electrode).trials(freq x time x trial)
        obs_pwr(:,:,:) = squeeze(ERSP{i_part,i_elect}.trials(:,timephi,:)); %get single subject's power
        % rank-transform power data (must be transformed)
        obs_pwr_reshape = reshape(obs_pwr,num_frex*nTimepoints,n_resp)';
        obs_pwr_rank = tiedrank(obs_pwr_reshape);
        

        % Technically, you want to perform a correlation, but a linear least-squares fit 
        % provides the same conceptual results while being ~15 times faster if you don't 
        % care about the actual scale of the data.
        % The following line shows how to compute the Spearman correlation coefficient:
        realcorrs_plot = 1-6*sum((obs_pwr_rank-repmat(errdeg_rank',1,size(obs_pwr_rank,2))).^2)/(n_resp*(n_resp^2-1));
        realcorrs_plot = reshape(realcorrs_plot,num_frex,nTimepoints);
        realcorrs_plot_out{i_part,ii} = realcorrs_plot; %save data
        
        % Convert rho values to t-vals
        tvalsRho{i_part,ii} = realcorrs_plot.*sqrt((n_resp-2)./(1-realcorrs_plot.^2));

        % The following line shows how to compute the Spearman correlation coefficient:
        %faster method for non-scaled data
        realcorrs = (errdeg_rank*errdeg_rank')\errdeg_rank*obs_pwr_rank; %faster
        realcorrs = reshape(realcorrs,num_frex,nTimepoints);
        realcorrs_out{i_part,ii} = realcorrs; %save data
        
        
        %% Create H0 distribution
        % initialize null hypothesis matrices
        permuted_rvals  = zeros(n_permutes,num_frex,nTimepoints);
        max_pixel_rvals = zeros(n_permutes,2);
        max_clust_info  = zeros(n_permutes,1);
        
        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:n_permutes
            fake_rt_mapping = errdeg_rank(randperm(n_resp));

            % compute non-scaled data correlations
            fakecorrs = (fake_rt_mapping*fake_rt_mapping')\fake_rt_mapping*obs_pwr_rank;

            % reshape to 2D map for cluster-correction
            fakecorrs = reshape(fakecorrs,num_frex,nTimepoints);

            % save all permuted values
            permuted_rvals(permi,:,:) = fakecorrs;

            % save maximum pixel values
            max_pixel_rvals(permi,:) = [ min(fakecorrs(:)) max(fakecorrs(:)) ];

            clear fakecorrs fake_rt_mapping
        end
        clear permi


        % now compute Z-map (standardized units from each correlation away from
        % the distribution of null-hypothesis correlation coefficients)
        zmap = (realcorrs-squeeze(mean(permuted_rvals,1)))./squeeze(std(permuted_rvals));
        zmap_out{i_part,ii} = zmap; %save values of plot
        
        
        %% two methods of evaluating statistical significance
        % Z-value to p %(p_z method)
        p_zcorr{i_part,ii} = 2*(1-normcdf(abs(zmap)));  %times 2 for 2-tailed test
%         p_zcorr{i_part,ii} = 2*normcdf(-abs(zmap));
        
        % p-value count
%         pCount{i_part,ii} = sum(permuted_rvals>realcorrs)/n_permutes;
        
        
        %% apply pixel-level corrected threshold
        lower_threshold = prctile(max_pixel_rvals(:,1),    mcc_voxel_pval*100/2);
        upper_threshold = prctile(max_pixel_rvals(:,2),100-mcc_voxel_pval*100/2);
        threshold_out{i_part,ii} = [lower_threshold upper_threshold];
        zmapthresh = zmap;
        zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=0;
        
        %p-values after correction (p_z method)
        p_zcorr_thresh{i_part,ii} = 2*(1-normcdf(abs(zmapthresh)));  %times 2 for 2-tailed test
        
        
        clear obs_pwr obs_pwr_reshape obs_pwr_rank realcorrs_plot permuted_rvals...
            max_pixel_rvals lower_threshold upper_threshold zmapthresh zmap...
            i_elect t
        
    end
    clear n_resp ii
    
end

clear errdeg errdeg_rank num_frex i_part nTimepoints voxel_pval...
    mcc_voxel_pval time_plot

% -------------------------------------------------------------------------
%% Save analyzed data
% some variables that are saved were not used in the paper
save([saveLocation 'corrR_bysubj_1e5perm.mat'],'freqs','p_tcorr','p_zcorr',...
    'p_zcorr_thresh','realcorrs','realcorrs_out','realcorrs_plot_out',...
    'threshold_out','times','tvalsRho','zmap_out','timewin'); 

% -------------------------------------------------------------------------

% clears variables that end/begin with...
clear -regexp \<n_


% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Plot of averages of z-scores

%extract data
thresh_tmp = NaN([length(exp.participants),length(exp.singletrialselecs),length(freqs),length(timephi)]);%pre-allocate
realcorrs_tmp = NaN([length(exp.participants),length(exp.singletrialselecs),length(freqs),length(times)]);%pre-allocate
realcorrs_plot_tmp = NaN([length(exp.participants),length(exp.singletrialselecs),length(freqs),length(timephi)]);%pre-allocate
tvalsRho_tmp = NaN([length(exp.participants),length(exp.singletrialselecs),length(freqs),length(timephi)]);%pre-allocate
zmap_tmp = NaN([length(exp.participants),length(exp.singletrialselecs),length(freqs),length(timephi)]);%pre-allocate
for ii = 1:length(exp.singletrialselecs)
    for i_part = 1:length(exp.participants) 
        thresh_tmp(i_part,ii,:,:) = threshold_out{i_part,ii};
        realcorrs_tmp(i_part,ii,:,:) = realcorrs_out{i_part,ii};
        realcorrs_plot_tmp(i_part,ii,:,:) = realcorrs_plot_out{i_part,ii};
        tvalsRho_tmp(i_part,ii,:,:) = tvalsRho{i_part,ii};
        zmap_tmp(i_part,ii,:,:) = zmap_out{i_part,ii};
    end
    clear i_part
end
clear ii i_elect

% Get mean of maps and threshold
cout_thresh = squeeze(mean(thresh_tmp,1));
cout_realcorrs = squeeze(mean(realcorrs_tmp,1));
cout_realcorrs_plot = squeeze(mean(realcorrs_plot_tmp,1));
cout_tvalsRho = squeeze(mean(tvalsRho_tmp,1));
cout_zmap = squeeze(mean(zmap_tmp,1));
clear thresh_tmp realcorrs_tmp realcorrs_plot_tmp tvalsRho_tmp zmap_tmp

% -------------------------------------------------------------------------
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
% v = df
% -------------------------------------------------------------------------

cmap = redblue(256); %create colormap colors
CLim = [-1 1]; %plot corr coefficients
tite = 'Corrs (Avg): ';
savname = 'Avg_zmap_';
for ii = 1:length(exp.singletrialselecs)

    corr_plot = squeeze(cout_zmap(ii,:,:));
    
    figure; colormap(cmap); %open a new figure
    imagesc(times(timephi),freqs,corr_plot,CLim)
    hold on
    
    % apply pixel-level corrected threshold
    zmapthresh = squeeze(cout_realcorrs(ii,:,:));
    lower_threshold = cout_thresh(ii,1);
    upper_threshold = cout_thresh(ii,2);
    zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=false;
    tmp_plot=logical(zmapthresh);
    
    % convert to p-values
%     tmp_plot = 2*(1-normcdf(abs(squeeze(cout_zmap(ii,:,:)))));  %times 2 for 2-tailed test
%     tmp_plot(tmp_plot>0.05)=false;
%     tmp_plot=logical(tmp_plot);

    % apply uncorrect threshold
%     zmapthresh = squeeze(cout_zmap(ii,:,:));
%     zmapthresh(abs(zmapthresh)<norminv(1-0.05))=0;
%     tmp_plot=logical(zmapthresh);
    
    contour(times,freqs,tmp_plot,1,'linecolor','k') %plot significant values
    
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    ylim([2 40]); yticks(5:5:40)
    xlim([-700 800]); xticks(-600:200:800)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Z-Score (Corrs)');
    t.Limits = CLim;
    title([tite exp.singtrlelec_name{ii}],'FontSize',10.5);
    
    savefig([saveFig savname exp.elec_names{ii}])
    
   clear i_elect i_cond tmp_plot zmapthresh lower_threshold upper_threshold corr_plot
end
clear ii ncond conds CLim cmap tite t savname


 
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Plot maps for each subject (not in manuscript)

cmap = redblue(256); %create colormap colors
CLim = [-.5 .5]; 

%choose data to plot
% plot_data = zmap_out; 
plot_data = realcorrs_plot_out;
% plot_data = tvalsRho;

%choose statistics
plot_cont = p_zcorr_thresh;
% plot_cont = p_zcorr;
% plot_cont = p_tcorr;

savname = 'pz_thresh_realcor_subj';

for i_part = 1:length(exp.participants) 
     
    figure('Position',[-1661 164 1492 747]); colormap(cmap); %open a new figure
    
    for ii = 1:length(exp.singletrialselecs)
        
        tmp_plot = plot_cont{i_part,ii}; %get data
        
        subtightplot(4,8,ii)
        imagesc(times(timephi),freqs,plot_data{i_part,ii},CLim) %plot z-map
        hold on
        
%         [h, crit_p, adj_ci_cvrg, tmp_plot] = fdr_bh(tmp_plot,0.05,'pdep','yes'); %multi-comp correction
        tmp_plot(tmp_plot>0.05)=false;
        tmp_plot=logical(tmp_plot);
        contour(times(timephi),freqs,tmp_plot,1,'linecolor','k','LineWidth',1)
        
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',0.75) %vertical line
        line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',0.75)  %vertical line for response screen onset
        ylim([2 40]); yticks(5:5:40)
        xlim([-700 800]); xticks(-600:200:800)
        ylabel('Freqency (Hz)'); xlabel('Time (ms)');
%         colorbar
%         t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
%         set(get(t,'ylabel'),'String', 'p-value');
        title(['sub ' num2str(i_part) ': ' exp.singtrlelec_name{ii}],'FontSize',9);
        
        savefig([saveFig savname num2str(i_part)])
        
        clear tmp_plot
    end
    clear n ii
end
 clear i_part plot_data cmap CLim plot_cont savname


 
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------






% clears variables that end/begin with...
clear -regexp \<n_ \<zmap_ _thresh\>











