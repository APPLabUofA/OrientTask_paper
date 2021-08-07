% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% This code uses data previously processed by Mdn_PwrR_Split.m and 
% Mdn_ERP_Split.m
% 
% The code below uses the model fitting function (fitlm) for the analysis 
% rather than the linear regression function (regress). This does not
% change the coefficients and R^2 statistics, but it does mean that the
% t-statistics, F-statistics and associated p-values are not accurate with
% regards to the analysis in the paper. To get the correct statistical
% test values, the regress function would also be needed along with the
% adjusted alpha level.
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
load('mdn_pwrR_split_v4.mat')

% Load data from Mdn_ERP_Split.m
load([saveLocation 'mdn_ERP_split_v4.mat'])

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;



% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% Guess Rate
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

% Combine into array
erpA_g(:,:,1) = g_aboveP1;
erpA_g(:,:,2) = g_aboveN1;
erpA_g(:,:,3) = g_aboveP2;
erpA_g(:,:,4) = g_aboveN2;
erpA_g(:,:,5) = g_aboveP3;

erpB_g(:,:,1) = g_belowP1;
erpB_g(:,:,2) = g_belowN1;
erpB_g(:,:,3) = g_belowP2;
erpB_g(:,:,4) = g_belowN2;
erpB_g(:,:,5) = g_belowP3;


% -------------------------------------------------------------------------
%% Get power separated data

clear timewin freqlim freqband

%finds the frequencies you want
freqband{3} = [8 40]; %alpha to gamma
freqband{2} = [4 7]; %theta
freqband{1} = [2 3]; %delta

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

% (participant x electrode x freq x ERP)
g_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqband),length(timewin)); %pre-allocate
g_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqband),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2)); %time range
            
            for ifreq = 1:length(freqband)
                freqlim = find(freqs>=(freqband{ifreq}(1)-0.5) & freqs<=(freqband{ifreq}(2)+0.5)); %freq range
                        % (participant x electrode x freq x time)
                g_Hpwr_stat(i_part,ii,ifreq,iwin) = squeeze(mean(mean(g_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
                g_Lpwr_stat(i_part,ii,ifreq,iwin) = squeeze(mean(mean(g_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            
                clear freqlim
            end
            clear ifreq
        end
        clear timelim iwin
    end
    clear ii i_elect
end
clear i_part timewin freqlim freqband


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Stepwise Linear Regression

% ERP_labels = {'P1';'N1';'P2';'N2';'P3'};
ERP_win    = [1,    2,   3,   4,   5]; %corresponds to the time windows

outH_g = cell(length(ERP_win),length(exp.singletrialselecs));
outL_g = cell(length(ERP_win),length(exp.singletrialselecs));
for ii = 1:length(ERP_win)
    
    for ielect = 1:length(exp.singletrialselecs)
    
        %//////////////////////////////////////////////////////////////////
        % High/Above
          % (participant x electrode x freq x ERP)
        dataX = squeeze(g_Hpwr_stat(:,ielect,1:3,ii));
        dataY = squeeze(erpA_g(:,ielect,ii));
   
        % Step 1
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        beta = (sx1/sy) * mdl.Coefficients.Estimate(2);
        clear sy sx1
        
        % Save results to data structure
        outH_g{ii,ielect}.b1 = mdl.Coefficients;
        outH_g{ii,ielect}.beta1 = beta;
        outH_g{ii,ielect}.r21 = mdl.Rsquared.Adjusted;
        clear beta mdl

        
        % Step 2
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        clear sy sx1 sx2
        
        % Save results to data structure
        outH_g{ii,ielect}.b2 = mdl.Coefficients;
        outH_g{ii,ielect}.beta2 = beta;
        outH_g{ii,ielect}.r22 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        
        % Step 3
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2 + x3');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        sx3 = std(mdl.Variables.x3);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        beta(3,1) = (sx3/sy) * mdl.Coefficients.Estimate(4);
        clear sy sx1 sx2 sx3
        
        % Save results to data structure
        outH_g{ii,ielect}.b3 = mdl.Coefficients;
        outH_g{ii,ielect}.beta3 = beta;
        outH_g{ii,ielect}.r23 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        clear dataX dataY
        
        
        %//////////////////////////////////////////////////////////////////
        % Low/Below
          % (participant x electrode x freq x ERP)
        dataX = squeeze(g_Lpwr_stat(:,ielect,1:3,ii));
        dataY = squeeze(erpB_g(:,ielect,ii));
        
        % Step 1
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        beta = (sx1/sy) * mdl.Coefficients.Estimate(2);
        clear sy sx1
        
        % Save results to data structure
        outL_g{ii,ielect}.b1 = mdl.Coefficients;
        outL_g{ii,ielect}.beta1 = beta;
        outL_g{ii,ielect}.r21 = mdl.Rsquared.Adjusted;
        clear beta mdl

        
        % Step 2
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        clear sy sx1 sx2
        
        % Save results to data structure
        outL_g{ii,ielect}.b2 = mdl.Coefficients;
        outL_g{ii,ielect}.beta2 = beta;
        outL_g{ii,ielect}.r22 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        
        % Step 3
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2 + x3');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        sx3 = std(mdl.Variables.x3);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        beta(3,1) = (sx3/sy) * mdl.Coefficients.Estimate(4);
        clear sy sx1 sx2 sx3
        
        % Save results to data structure
        outL_g{ii,ielect}.b3 = mdl.Coefficients;
        outL_g{ii,ielect}.beta3 = beta;
        outL_g{ii,ielect}.r23 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        
        clear dataX dataY
  
    end
    clear ielect

end
clear ii

% -------------------------------------------------------------------------
% Save variables to .mat file
save([saveLocation 'regress_ERP_PwrR.mat'],'ERP_labels','erpA_g','erpB_g',...
    'g_Hpwr_stat','g_Lpwr_stat','outH_g','outL_g')
% -------------------------------------------------------------------------

clear nperm tail alpha stat timewin freqlim freqband




% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
%% Standard Deviation Parameter
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------

% Combine into array
erpA_sd(:,:,1) = sd_aboveP1;
erpA_sd(:,:,2) = sd_aboveN1;
erpA_sd(:,:,3) = sd_aboveP2;
erpA_sd(:,:,4) = sd_aboveN2;
erpA_sd(:,:,5) = sd_aboveP3;

erpB_sd(:,:,1) = sd_belowP1;
erpB_sd(:,:,2) = sd_belowN1;
erpB_sd(:,:,3) = sd_belowP2;
erpB_sd(:,:,4) = sd_belowN2;
erpB_sd(:,:,5) = sd_belowP3;


% -------------------------------------------------------------------------
%% Get power separated data

clear timewin freqlim freqband

%finds the frequencies you want
freqband{3} = [8 40]; %alpha to gamma
freqband{2} = [4 7]; %theta
freqband{1} = [2 3]; %delta

%finds the times you want from the timess variable
timewin{1} = [80 140]; %P1
timewin{2} = [140 200]; %N1
timewin{3} = [200 255]; %P2
timewin{4} = [255 360]; %N2
timewin{5} = [360 500]; %P3

% (participant x electrode x freq x ERP)
sd_Hpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqband),length(timewin)); %pre-allocate
sd_Lpwr_stat = NaN(length(exp.participants),length(exp.singletrialselecs),length(freqband),length(timewin)); %pre-allocate
for i_part = 1:length(exp.participants)
    for ii = 1:length(exp.singletrialselecs)
        i_elect = exp.singletrialselecs(ii); %for doing only a selection of electrodes
        
        for iwin = 1:length(timewin)
            timelim = find(times>=timewin{iwin}(1) & times<=timewin{iwin}(2)); %time range
            
            for ifreq = 1:length(freqband)
                freqlim = find(freqs>=(freqband{ifreq}(1)-0.5) & freqs<=(freqband{ifreq}(2)+0.5)); %freq range
                        % (participant x electrode x freq x time)
                sd_Hpwr_stat(i_part,ii,ifreq,iwin) = squeeze(mean(mean(sd_out_Hpwr(i_part,i_elect,freqlim,timelim),4),3)); 
                sd_Lpwr_stat(i_part,ii,ifreq,iwin) = squeeze(mean(mean(sd_out_Lpwr(i_part,i_elect,freqlim,timelim),4),3));
            
                clear freqlim
            end
            clear ifreq
        end
        clear timelim iwin
    end
    clear ii i_elect
end
clear i_part timewin freqlim freqband


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Stepwise Linear Regression

% ERP_labels = {'P1';'N1';'P2';'N2';'P3'};
ERP_win    = [1,    2,   3,   4,   5]; %corresponds to the time windows

outH_sd = cell(length(ERP_win),length(exp.singletrialselecs));
outL_sd = cell(length(ERP_win),length(exp.singletrialselecs));
for ii = 1:length(ERP_win)
    
    for ielect = 1:length(exp.singletrialselecs)
        
        %//////////////////////////////////////////////////////////////////
        % High/Above
          % (participant x electrode x freq x ERP)
        dataX = squeeze(sd_Hpwr_stat(:,ielect,1:3,ii));
        dataY = squeeze(erpA_sd(:,ielect,ii));
   
        % Step 1
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        beta = (sx1/sy) * mdl.Coefficients.Estimate(2);
        clear sy sx1
        
        % Save results to data structure
        outH_sd{ii,ielect}.b1 = mdl.Coefficients;
        outH_sd{ii,ielect}.beta1 = beta;
        outH_sd{ii,ielect}.r21 = mdl.Rsquared.Adjusted;
        clear beta mdl

        
        % Step 2
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        clear sy sx1 sx2
        
        % Save results to data structure
        outH_sd{ii,ielect}.b2 = mdl.Coefficients;
        outH_sd{ii,ielect}.beta2 = beta;
        outH_sd{ii,ielect}.r22 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        
        % Step 3
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2 + x3');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        sx3 = std(mdl.Variables.x3);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        beta(3,1) = (sx3/sy) * mdl.Coefficients.Estimate(4);
        clear sy sx1 sx2 sx3
        
        % Save results to data structure
        outH_sd{ii,ielect}.b3 = mdl.Coefficients;
        outH_sd{ii,ielect}.beta3 = beta;
        outH_sd{ii,ielect}.r23 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        clear dataX dataY
        
        
        %//////////////////////////////////////////////////////////////////
        % Low/Below
          % (participant x electrode x freq x ERP)
        dataX = squeeze(sd_Lpwr_stat(:,ielect,1:3,ii));
        dataY = squeeze(erpB_sd(:,ielect,ii));
        
        % Step 1
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        beta = (sx1/sy) * mdl.Coefficients.Estimate(2);
        clear sy sx1
        
        % Save results to data structure
        outL_sd{ii,ielect}.b1 = mdl.Coefficients;
        outL_sd{ii,ielect}.beta1 = beta;
        outL_sd{ii,ielect}.r21 = mdl.Rsquared.Adjusted;
        clear beta mdl

        
        % Step 2
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        clear sy sx1 sx2
        
        % Save results to data structure
        outL_sd{ii,ielect}.b2 = mdl.Coefficients;
        outL_sd{ii,ielect}.beta2 = beta;
        outL_sd{ii,ielect}.r22 = mdl.Rsquared.Adjusted;
        clear beta mdl
        
        
        % Step 3
        mdl = fitlm(dataX,dataY,'y ~ 1 + x1 + x2 + x3');

        % Calculate beta weights/standardized beta
        sy = std(mdl.Variables.y);
        sx1 = std(mdl.Variables.x1);
        sx2 = std(mdl.Variables.x2);
        sx3 = std(mdl.Variables.x3);
        beta(1,1) = (sx1/sy) * mdl.Coefficients.Estimate(2);
        beta(2,1) = (sx2/sy) * mdl.Coefficients.Estimate(3);
        beta(3,1) = (sx3/sy) * mdl.Coefficients.Estimate(4);
        clear sy sx1 sx2 sx3
        
        % Save results to data structure
        outL_sd{ii,ielect}.b3 = mdl.Coefficients;
        outL_sd{ii,ielect}.beta3 = beta;
        outL_sd{ii,ielect}.r23 = mdl.Rsquared.Adjusted;
        clear beta mdl
        

        clear dataX dataY
    end
    clear ielect

end
clear ii

% -------------------------------------------------------------------------
% Save variables to existing .mat file
save([saveLocation 'regress_ERP_PwrR.mat'],'erpA_sd','erpB_sd',...
    'sd_Hpwr_stat','sd_Lpwr_stat','outH_sd','outL_sd','-append')
% -------------------------------------------------------------------------

clear nperm tail alpha stat timewin freqlim freqband ERP_win ERP_labels





% /////////////////////////////////////////////////////////////////////////
% Clear workspace
ccc
% /////////////////////////////////////////////////////////////////////////













