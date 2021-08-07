
% /////////////////////////////////////////////////////////////////////////
%% Load settings
load('filt_byTargets_v4_Settings.mat'); %setting for target-aligned trials

% /////////////////////////////////////////////////////////////////////////
%% Location of saved ERP results
saveLocation = [exp.dataLocation '\ProcessData\']; % set save directory of data set

% /////////////////////////////////////////////////////////////////////////
%% Location to save figures
saveFig = [pwd '\Figures\ERP\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end


% /////////////////////////////////////////////////////////////////////////
%% ------------------ Load Previously Saved Data --------------------------
% /////////////////////////////////////////////////////////////////////////

% Load data from ERP_byTarget.m
load([saveLocation 'erp_out_byTarget.mat'])

% Load data from Anal_ERP_Orient_byCatchTrial.m
erp_catchtrials = cell2mat(struct2cell(load([saveLocation 'erp_out_catch_trials.mat'],'erp_out')));

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Select electrodes for analysis

% List electrodes with EOGs if want to check EOGs
% elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2';'VEOG';'HEOG'};

% List electrodes to get ERPs and topograph plots 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% ::::::::::::::::::  Plot the ERPs by electrode  ::::::::::::::::::::::::
% :::::::::::::::::::     Subtract Catch Trial     ::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% ------------------------------------------------------------------------- 
% Average across subjects by errors
erp_out_byerr2(:,:,1) = squeeze((mean(erp_out_x(:,:,:),1)-mean(erp_catchtrials(:,:,:),1))); %small errors
erp_out_byerr2(:,:,2) = squeeze((mean(erp_out_n(:,:,:),1)-mean(erp_catchtrials(:,:,:),1))); %large errors
% Average across subjects - catch trials
erp_out_bycatch(:,:) = squeeze(mean(erp_catchtrials(:,:,:),1));
% ------------------------------------------------------------------------- 

% Plot ERPs with error bars
for ii = 1:length(elect_erp) 
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
%     i_elect = ii; %if want to select one making ERPs

    % get axes limits
    xmin = -200; xmax = 800;
    ymin = -6; ymax = 6; %brain electrodes
%     ymin = -5; ymax = 30; %VEOG: ii=33
%     ymin = -2; ymax = 4; %HEOG: ii=34
    
    figure('Color',[1 1 1]); 
    boundedline(erp_time,erp_out_byerr2(i_elect,:,1),...
        squeeze(std((erp_out_x(:,i_elect,:)-erp_catchtrials(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'c',...
            erp_time,erp_out_byerr2(i_elect,:,2),...
            squeeze(std((erp_out_n(:,i_elect,:)-erp_catchtrials(:,i_elect,:)),[],1))./sqrt(length(exp.participants)),'m');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    
    title([el_erp_names{ii} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); 
    xticks(-200:100:xmax); 
    yticks(ymin:2:ymax)
%     grid on
        
    legend({'Accurate','Guess'},'Location','best');
    
    savefig([saveFig 'ERP-NoT_' el_erp_names{ii}])
end
clear ii xmax xmin ymin ymax   

% ------------------------------------------------------------------------- 



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


% List electrodes to get topograph plots 
elect_erp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% el_erp_names = {'M2';'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
%     'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% Difference topplots ---------------------------------------------------------------------------- 
% average across subjects by errors
erp_out_byerr(:,:,1) = squeeze(mean((erp_out_x(:,:,:)-erp_catchtrials(:,:,:)),1)); %small errors
erp_out_byerr(:,:,2) = squeeze(mean((erp_out_n(:,:,:)-erp_catchtrials(:,:,:)),1)); %large errors
erp_out_diff(:,:,1) = erp_out_byerr(:,:,1) - erp_out_byerr(:,:,2); %small-large
% average across subjects - catch trials
erp_out_bycatch(:,:) = squeeze(mean(erp_catchtrials(:,:,:),1));

clear erp_out_byerr
% ------------------------------------------------------------------------------------------------ 

% Set the range of time to consider
tWin{1} = [80 140]; %P1
tWin{2} = [140 200]; %N1
tWin{3} = [200 255]; %P2
tWin{4} = [255 360]; %N2
tWin{5} = [360 500]; %P3

CLims1 = [-3 3]; %range in microvolts
nconds = 1; %number of plots
for tw_i = 1:length(tWin) %loop through several time windows 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %this code finds the times you want from the times variable
    time_window = find(erp_time>= itWin(1),1):find(erp_time>= itWin(2),1)-1;
    
    figure('Color',[1 1 1]); %set-up figure
    
    %loop through conditions to make plot of each
    for i_cond = 1:nconds         
        set(gca,'Color',[1 1 1]);
        
        temp = mean(erp_out_diff(:,time_window,i_cond),2)'; %ERP within time window
        temp(1) = NaN; %so M2 is not included

        topoplot(temp,chan_locs,'whitebk','on','plotrad',0.6,'maplimits',CLims1,...
        'plotchans',elect_erp,'emarker',{'.','k',11,1})
        title('Accurate-Guess');
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage (uV)');
        clear temp
    end
    % Overall subplot title
    supertitle([num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)
    
    savefig([saveFig 'ERP-NoT_Topo_' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'])
   
    clear itWin time_window i_cond

end
clear tw_i nconds conds tWin CLims CLims1 t time1 time2 i_elect


clear erp_out_diff



% ------------------------------------------------------------------------- 
% /////////////////////////////////////////////////////////////////////////
%% :::::::::::::::::::::  Permutation Test  :::::::::::::::::::::::::::::::
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ------------------------------------------------------------------------- 
% **Need code from: https://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox

% ------------------------------------------------------------------------- 
% Just brain electrodes
elect_erp = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
el_erp_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

% ------------------------------------------------------------------------- 

% ERPs in paper
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

% electrode x time points x conditions x subjects array 
erp_stat = NaN(length(elect_erp),length(timewin),2,length(exp.participants)); %pre-allocate
for i_time = 1:length(timewin)
    timelim = find(erp_time>=timewin{i_time}(1) & erp_time<=timewin{i_time}(2));

    for i_part = 1:length(exp.participants)
        for ii = 1:length(elect_erp)
            i_elect = elect_erp(ii); %for doing only a selection of electrodes        
                                                 % (participant x electrode x time)
            % Subtract out catch trials
            erp_stat(i_elect,i_time,1,i_part) = squeeze(mean((erp_out_x(i_part,i_elect,timelim)-erp_catchtrials(i_part,i_elect,timelim)),3)); %small errors
            erp_stat(i_elect,i_time,2,i_part) = squeeze(mean((erp_out_n(i_part,i_elect,timelim)-erp_catchtrials(i_part,i_elect,timelim)),3)); %large errors
           
            % Save electrode names
            permtest.elect_id{ii,1} = el_erp_names{ii};
        end
        clear ii i_elect
    end
    clear i_part timelim
end
clear i_time

% One sample/repeated-measures permutation test
% Open pval for outcome 
% data(Channel x Time x Participant)
[pval, t_orig, tmx_ptile] = mxt_perm1(squeeze(erp_stat(2:32,:,1,:)-erp_stat(2:32,:,2,:)),1e5);

writematrix(pval,'ERP_permtest_pval.txt','Delimiter','tab')

clear pval t_orig tmx_ptile seed_state est_alpha timewin


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% clears variables that end/begin with...
clear -regexp \<permtest_ \<ttest_
clear time_window time1 time2 erp_stat 

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%% ::::::::::::::::::  Plot the ERPs by electrode  ::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ------------------------------------------------------------------------- 
% Average across subjects by errors
erp_out_byerr(:,:,1) = squeeze(mean(erp_out_x(:,:,:),1)); %small errors
erp_out_byerr(:,:,2) = squeeze(mean(erp_out_n(:,:,:),1)); %large errors
erp_out_byerr(:,:,3) = squeeze(mean(erp_catchtrials(:,:,:),1)); %catch trials (no target)
% Average across subjects - catch trials
erp_out_bycatch(:,:) = squeeze(mean(erp_catchtrials(:,:,:),1));
% ------------------------------------------------------------------------- 


% Plot ERPs with error bars
for ii = 1:length(elect_erp) 
    i_elect = elect_erp(ii); %for doing only a selection of electrodes
%     i_elect = ii; %selection done when making ERPs
    % get axes limits
    xmin = -200; xmax = 800;
    ymin = -6; ymax = 10;
%     ymin = -5; ymax = 30; %VEOG: ii=33
%     ymin = -2; ymax = 4; %HEOG: ii=34
    
    figure('Color',[1 1 1]); 
    boundedline(erp_time,erp_out_bycatch(i_elect,:),squeeze(std(erp_catchtrials(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'k',...
            erp_time,erp_out_byerr(i_elect,:,1),squeeze(std(erp_out_x(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'c',...
            erp_time,erp_out_byerr(i_elect,:,2),squeeze(std(erp_out_n(:,i_elect,:),[],1))./sqrt(length(exp.participants)),'m');
    hold on
    line([xmin xmax],[0 0],'color','k','LineWidth',1.5) %horizontal line
    line([0 0],[ymin ymax],'color','k','LineWidth',1.5) %vertical line
%     line([50 50],[ymin ymax],'LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[ymin ymax],'color','r','LineStyle','--','LineWidth',1.5)  %vertical line for color wheel onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    
    title([el_erp_names{ii} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); 
    xticks(-200:100:xmax); 
    yticks(ymin:2:ymax)
%     grid on
        
    legend({'No Target','Accurate','Guess'},'Location','best');
    
    savefig([saveFig 'ERP_' el_erp_names{ii}])

end
clear ii xmax xmin ymin ymax   

% ------------------------------------------------------------------------- 





