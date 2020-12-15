% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% Does multiple runs of the nonparametric permutation version of Spearman's
% correlation between fitted parameter values and power across the
% time-frequency space. Uses pixel-based multiple corrections procedure.
% Saves the runs to create an average of the results to ensure stability.
% Use Permut_Corr_AvgMaps.m to combine results and make plots
% 
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Location to save permutation results
saveLocation_perm = [exp.dataLocation '\ProcessData\permut_data\']; 
saveFile_perm = 'Corr_SD_Pwr\';

% if folder doesn't exist yet, create one
if ~exist([saveLocation_perm saveFile_perm])
    mkdir([saveLocation_perm saveFile_perm]);
end

%% Load previously processed target-aligned epoch data
%data created in Power_ZScore_ERS_Plots.m
all_ersp_Z = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'all_ersp_z'));  %gets loaded as a struct
all_ersp_Z = all_ersp_Z{1};
chanlocs = struct2cell(load([saveLocation 'all_ersp_Z_v4.mat'],'chanlocs')); 
chanlocs = chanlocs{1};

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Extract standard mixed model parameter values
% loaded from saved behavioral data or can be created in BEH_ModelFits_OrientTask.m

sd_out_cat = NaN(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants)
%     g_out_cat(i_part) = model_out{1,i_part}(1);
    sd_out_cat(i_part) = model_out{1,i_part}(2);
end
clear i_part

% rank-transform parameter
SD_rank = tiedrank(sd_out_cat);

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

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Permutation Testing
%set parameters inside for loop
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction

n_permutes = 1e5; %number of permutations used to estimate null distribution

num_frex = length(freqs); %number of frequencies
nTimepoints = length(times); %number of timepoints
n_resp = length(SD_rank); %total number of subjects

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

%for repeating analysis & saving results to make an avg map
for reps = 1:10

    zmap_out = cell(1,length(exp.singletrialselecs)); %pre-allocate
    threshold_out = cell(1,length(exp.singletrialselecs)); %pre-allocate
    realcorrs_out = cell(1,length(exp.singletrialselecs)); %pre-allocate
    realcorrs_plot_out = cell(1,length(exp.singletrialselecs)); %pre-allocate

    for ii = 1:length(exp.singletrialselecs)

        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes

        % rank-transform power data (must be transformed)
        obs_pwr_reshape = reshape(obs_pwr{i_elect},num_frex*nTimepoints,n_resp)';
        obs_pwr_rank = tiedrank(obs_pwr_reshape);

        % Technically, you want to perform a correlation, but a linear least-squares fit 
        % provides the same conceptual results while being ~15 times faster if you don't 
        % care about the actual scale of the data.
        % The following line shows how to compute the Spearman correlation coefficient:
        realcorrs_plot = 1-6*sum((obs_pwr_rank-repmat(SD_rank',1,size(obs_pwr_rank,2))).^2)/(n_resp*(n_resp^2-1));
        realcorrs_plot = reshape(realcorrs_plot,num_frex,nTimepoints);
        realcorrs_plot_out{ii} = realcorrs_plot; %save data
        %faster method for non-scaled data
        realcorrs = (SD_rank*SD_rank')\SD_rank*obs_pwr_rank; %faster
        realcorrs = reshape(realcorrs,num_frex,nTimepoints);
        realcorrs_out{ii} = realcorrs; %save data

        % initialize null hypothesis matrices
        permuted_rvals  = zeros(n_permutes,num_frex,nTimepoints);
        max_pixel_rvals = zeros(n_permutes,2);
        max_clust_info  = zeros(n_permutes,1);

        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:n_permutes
            fake_rt_mapping = SD_rank(randperm(n_resp));

            % compute t-map of null hypothesis
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
        zmap_out{ii} = zmap; %save values of plot

        % apply pixel-level corrected threshold
        lower_threshold = prctile(max_pixel_rvals(:,1),    mcc_voxel_pval*100/2);
        upper_threshold = prctile(max_pixel_rvals(:,2),100-mcc_voxel_pval*100/2);
        threshold_out{ii} = [lower_threshold upper_threshold];

        zmapthresh = zmap;
        zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=false;
        zmapthresh=logical(zmapthresh);

        clear zmap zmapthresh realcorrs permuted_rvals max_pixel_rvals...
            max_clust_info whichclusters2remove clust_info clustinfo clust_threshold...
            lower_threshold upper_threshold obs_pwr_reshape obs_pwr_rank realcorrs_plot

    end
    clear ii i_elect n t real_condition_mapping

    % Save each data file
    save_file = [saveLocation_perm saveFile_perm 'zmap_corr_SD_' num2str(reps) '.mat'];
    save(save_file,'threshold_out','zmap_out','realcorrs_out','realcorrs_plot_out')
    
    clear threshold_out zmap_out realcorrs_out save_file realcorrs_plot_out

end


% -------------------------------------------------------------------------

% clears variables that end/begin with...
clear -regexp \<n_ \<zmap_ _thresh\>

clear errdeg errdeg_rank reps num_frex SD_rank sd_out_cat...
    nTimepoints mcc_voxel_pval saveLocation_perm saveFile_perm

% -------------------------------------------------------------------------























