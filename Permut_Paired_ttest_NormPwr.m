% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Location to save permutation results
saveLocation_perm = [exp.dataLocation '\ProcessData\permut_data\']; 
saveFile_perm = 'paired_ttest_pwr\';

% if folder doesn't exist yet, create one
if ~exist([saveLocation_perm saveFile_perm])
    mkdir([saveLocation_perm saveFile_perm]);
end

%% Load previously processed target-aligned epoch data
%data comes from Power_ZScore_ERS_Plots.m
load([saveLocation 'pwr_AvG_v4.mat']);

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Permutation Testing

%finds the times you want from the times variable
timewin = [-700 800];
timephi = find(times>=timewin(1) & times<=timewin(2));


%set test parameters inside for loop
for reps = 1:10 %for repeating analysis & saving results to make an avg map

    voxel_pval = 0.05;
    mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
    n_permutes = 1e5; %number of permutations

    n = 24; %number of subjects
    num_frex = length(freqs); %number of frequencies

    nTimepoints = length(timephi); %number of timepoints
    
    % if folder doesn't exist yet, create one
    if ~exist([saveLocation_perm saveFile_perm])
        mkdir([saveLocation_perm saveFile_perm]);
    end

    zmap_out = cell(1,length(exp.singletrialselecs)); %pre-allocate
    threshold_out = cell(1,length(exp.singletrialselecs)); %pre-allocate
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.electrode(ii); %get electrode number

        % Put power data into one variable
        % (participants x electrodes x frequencies x timepoints)
        tmp_eegpwr = squeeze(cat(1,x_pwr{i_elect}(:,:,timephi),n_pwr{i_elect}(:,:,timephi)));

        % Logical to select data for t-test computation
        real_condition_mapping = [-ones(1,n) ones(1,n)];

        % compute actual paired t-test of difference
        % note. one-sample test of difference is = to paired t-test
        tnum   = squeeze(nanmean(tmp_eegpwr(real_condition_mapping==-1,:,:),1) - nanmean(tmp_eegpwr(real_condition_mapping==1,:,:),1));
        tdenom = squeeze(std((tmp_eegpwr(real_condition_mapping==-1,:,:) - tmp_eegpwr(real_condition_mapping==1,:,:)),[],1) ./ sqrt(n));
        real_t = tnum./tdenom;
        clear tnum tdenom

        % initialize null hypothesis matrices
        permuted_tvals  = zeros(n_permutes,num_frex,nTimepoints);
        max_pixel_pvals = zeros(n_permutes,2);
        max_clust_info  = zeros(n_permutes,1);

        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:n_permutes
        %     fake_condition_mapping = sign(randn(n*2,1));
            fake_condition_mapping = real_condition_mapping(randperm(n*2)); %need equal number of data points in both groups

            % compute t-map of null hypothesis
            tnum   = squeeze(nanmean(tmp_eegpwr(fake_condition_mapping==-1,:,:),1) - nanmean(tmp_eegpwr(fake_condition_mapping==1,:,:),1));
            tdenom = squeeze(std((tmp_eegpwr(fake_condition_mapping==-1,:,:) - tmp_eegpwr(fake_condition_mapping==1,:,:)),[],1) ./ sqrt(n));
            tmap   = tnum./tdenom;

            % save all permuted values
            permuted_tvals(permi,:,:) = tmap;

            % save maximum pixel values
            max_pixel_pvals(permi,:) = [ min(tmap(:)) max(tmap(:)) ];

            % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
            % note that here, clusters were obtained by parametrically thresholding
            % the t-maps
            tmap(abs(tmap)<tinv(1-voxel_pval,n-1))=0;

            % get number of elements in largest supra-threshold cluster
            clustinfo = bwconncomp(tmap);
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero

            clear tnum tdenom fake_condition_mapping
        end
        clear permi

        % now compute Z-map
        zmap = (real_t-squeeze(nanmean(permuted_tvals,1)))./squeeze(std(permuted_tvals));
        zmap_out{i_elect} = zmap; %save values of plot


        % apply pixel-level corrected threshold
        lower_threshold = prctile(max_pixel_pvals(:,1),    mcc_voxel_pval*100/2);
        upper_threshold = prctile(max_pixel_pvals(:,2),100-mcc_voxel_pval*100/2);
        threshold_out{i_elect} = [lower_threshold upper_threshold];

        zmapthresh = zmap;
        zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=0;


        clear tmp_eegpwr lower_threshold upper_threshold zmap zmapthresh tmap permuted_tvals...
            permuted_tvals max_pixel_pvals max_clust_info real_t whichclusters2remove...
            clust_info clustinfo clust_threshold
    end
    clear ii n_permutes i_elect n t num_frex nTimepoints voxel_pval mcc_cluster_pval...
        mcc_voxel_pval cmap real_condition_mapping CLim
    
    
    % Save each data file
    save_file = [saveLocation_perm saveFile_perm 'zmap_pwr_' num2str(reps) '.mat'];
    save(save_file,'threshold_out','zmap_out','timewin')


    clear threshold_out zmap_out save_file

end

clear reps saveLocation_perm saveFile_perm timephi timewin






