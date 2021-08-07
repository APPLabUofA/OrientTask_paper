% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                               INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% Need the MemToolbox found here: http://visionlab.github.io/MemToolbox/
% Fitting to mixed model and goodness-of-fit are done with functions from
% the MemToolbox

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load settings
load('filt_byTargets_v4_Settings.mat');

%% General location of saved processed data
saveLocation = [exp.dataLocation '\ProcessData\']; 

%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% OR load data file with response errors for each subject (and other task data)
%this loads all behavioral data including trials that might be excluded due
%to bad EEG data
error_deg = cell(1,length(exp.participants)); %pre-allocate
for i_part = 1:length(exp.participants)
    
    % Load each subject's data
    load([exp.dataLocation '\RawData\BEH\' num2str(exp.participants{i_part}) '_Orient_v3.mat'])

    % Remove practice trials (first 20 trials)
    error_deg_tmp = data.errorDegrees(1,21:end);
    targets = data.target(1,21:length(data.errorDegrees));

    % Remove trials where no target appeared
    error_deg{i_part} = error_deg_tmp(targets == 1);
    
    clear error_deg_tmp targets data
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Compare fits of errors to two different model
% Akaike Information Criterion (AIC): + values evidence in favor of model 2, - in favor of model 1
% **Bayesian Information Criterion (BIC): + values evidence in favor of model 2, - in favor of model 1
% Log likelihood: - values evidence in favor of model 2, + in favor of model 1 
% Corrected Akaike Information (AICc): + values evidence in favor of model 2, - in favor of model 1
% Deviance Information Criterion (DIC): + values evidence in favor of model 2, - in favor of model 1

model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_cmp = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
    sprintf(['Subj ' num2str(ii)])
    model_cmp{ii} = MemFit(resp_errdeg{ii},{model1,model2}); %comparing models
    
    %Make differences for easy comparison
    cmpDiff.AIC(ii,1) = model_cmp{ii}.AIC(1,1) - model_cmp{ii}.AIC(1,2);
    cmpDiff.BIC(ii,1) = model_cmp{ii}.BIC(1,1) - model_cmp{ii}.BIC(1,2);
    cmpDiff.logLike(ii,1) = model_cmp{ii}.logLike(1,1) - model_cmp{ii}.logLike(1,2);
    cmpDiff.AICc(ii,1) = model_cmp{ii}.AICc(1,1) - model_cmp{ii}.AICc(1,2);
end
clear ii error_deg model1 model2
T = struct2table(cmpDiff);

% Plot fit differences
figure; boxplot([T.AIC,T.BIC,T.AICc],'Labels',{'Akaike Information Criterion (AIC)',...
    'Bayesian Information Criterion (BIC)','Corrected Akaike Information (AICc)'},'Whisker',2.5)
ylim([-10 30])
ylabel(sprintf('Model Comparison Metric \n (Standard Mixture vs Standard Mixture w/Bias)'));


clear cmpDiff


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Fit errors to standard mixed model with plotting
model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_out = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
    if ii == 21 || ii == 22 %subjects had bias in mean errors
        model_out{ii} = MemFit(resp_errdeg{ii},model2); %plots fits  
%         model_out{ii} = MLE(resp_errdeg{ii},model2); %fits without plotting - w/bias
        model_out{ii} = model_out{ii}(2:3);
    else
        model_out{ii} = MemFit(resp_errdeg{ii},model1); %plots fits
%         model_out{ii} = MLE(resp_errdeg{ii},model1); %fits without plotting
    end
end
clear ii error_deg model model1 model2

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Figure 2A in manuscript
% Very annoying to do...
% 1) Select one participants plot at a time
% 2) Run function below
        fig=gco;
% 3) Save each participants' y-values in variable below. Change number in
% dataY_All(#,:) to correspond to the participant
        dataY_All(26,:)=fig.YData;
% 4) Save the x-values once since they are the same for each participant
        dataX_All(1,:)=fig.XData;

        
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Definitely want to save data so do not have to do again        
save([saveLocation 'BEH_Fits_All_v4.mat'],'dataY_All','dataX_All')

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Boundline plot
figure('Color',[1 1 1]); 
boundedline(dataX_All(1,:),squeeze(mean(dataY_All(:,:),1)),...
    squeeze(std(dataY_All(:,:),[],1)./sqrt(26)),'k')
xlim([-180 180])
title(['Overall Model Fits (Mean +/- SEM)']);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Put all the response errors across subjects into vector
resp_errdeg_cat = cat(2,resp_errdeg{1:end});

% Fit all errors to mixed model
model = StandardMixtureModel(); %standard 2 parameter model
% model = WithBias(StandardMixtureModel); %model with mu
model_out_cat = MemFit(resp_errdeg_cat,model); %fits with plotting
clear model

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Figure 2B in manuscript
% Extract parameter values
g_out_cat = NaN(length(exp.participants),1);
sd_out_cat = NaN(length(exp.participants),1);
for i_part = 1:length(exp.participants)
    g_out_cat(i_part) = model_out{1,i_part}(1);
    sd_out_cat(i_part) = model_out{1,i_part}(2);
end
clear i_part

median(g_out_cat)
std(g_out_cat)
mean(g_out_cat)

median(sd_out_cat)
mean(sd_out_cat)
std(sd_out_cat)

% Boxplots
figure
boxplot(g_out_cat,'Labels',{'g'},'Whisker',2)
ylabel('Guess Rate'); 
ylim([0 0.8])

figure
boxplot(sd_out_cat,'Labels',{'SD'},'Whisker',2)
ylabel('Variability'); 
ylim([5 20])


clear g_out_cat sd_out_cat


% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////








