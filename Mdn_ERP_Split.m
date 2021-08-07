% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% This code uses data previously processed by Mdn_PwrR_Split.m
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
% Load data from ERP_byTarget.m
load([saveLocation 'erp_out_all_v4.mat'])

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data and time & freq parameters
load([exp.dataLocation '\ProcessData\ALLEEG_' exp.settings '.mat'])

%initialize EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;




% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% P1 ERP
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

time_win = [80 140]; %P1
timelim = find(times_erp>time_win(1),1)-1:find(times_erp>time_win(2),1)-1;

errdeg_aboveP1 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
errdeg_belowP1 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
g_aboveP1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
g_belowP1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_aboveP1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_belowP1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for selecting specific electrodes
        
        erp_tmp = squeeze(mean(erp_out_all{1,i_part}(i_elect,timelim,:),2)); %average over time
        cutoff = median(erp_tmp); %get median to split trials
        
        above_mdn = find(erp_tmp > cutoff);
        below_mdn = find(erp_tmp <= cutoff);
        
        errdeg_above = resp_errdeg{i_part}(above_mdn);
        errdeg_below = resp_errdeg{i_part}(below_mdn);
        
        if i_part == 21 || i_part == 22 %subjects had bias in mean errors
            model_out_tmp = MLE(errdeg_above,model2); %fits without plotting
            model_out_above = model_out_tmp(2:3);
            clear model_out_tmp
            
            model_out_tmp = MLE(errdeg_below,model2); %fits without plotting
            model_out_below = model_out_tmp(2:3);
            clear model_out_tmp
        else
            model_out_above = MLE(errdeg_above,model); %fits without plotting
            model_out_below = MLE(errdeg_below,model); %fits without plotting
        end

        g_aboveP1(i_part,ii) = model_out_above(1);
        g_belowP1(i_part,ii) = model_out_below(1);
        sd_aboveP1(i_part,ii) = model_out_above(2);
        sd_belowP1(i_part,ii) = model_out_below(2);
        
        errdeg_aboveP1{i_part,ii} = errdeg_above;
        errdeg_belowP1{i_part,ii} = errdeg_below;
        
        clear model_out_above model_out_below cutoff above_mdn below_mdn...
            erp_tmp i_elect errdeg_above errdeg_below
    end
    clear ii
end
clear i_part model model2 time_win timelim


% Save variables to .mat file
save([saveLocation 'mdn_ERP_split_v4.mat'],'errdeg_aboveP1','errdeg_belowP1','g_aboveP1',...
    'g_belowP1','sd_aboveP1','sd_belowP1','times_erp','chan_locs')


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% N1 ERP
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

time_win = [140 200]; %N1
timelim = find(times_erp>time_win(1),1)-1:find(times_erp>time_win(2),1)-1;

errdeg_aboveN1 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
errdeg_belowN1 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
g_aboveN1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
g_belowN1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_aboveN1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_belowN1 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for selecting specific electrodes
        
        erp_tmp = squeeze(mean(erp_out_all{1,i_part}(i_elect,timelim,:),2)); %average over time
        cutoff = median(erp_tmp); %get median to split trials
        
        above_mdn = find(erp_tmp > cutoff);
        below_mdn = find(erp_tmp <= cutoff);
        
        errdeg_above = resp_errdeg{i_part}(above_mdn);
        errdeg_below = resp_errdeg{i_part}(below_mdn);
        
        if i_part == 21 || i_part == 22 %subjects had bias in mean errors
            model_out_tmp = MLE(errdeg_above,model2); %fits without plotting
            model_out_above = model_out_tmp(2:3);
            clear model_out_tmp
            
            model_out_tmp = MLE(errdeg_below,model2); %fits without plotting
            model_out_below = model_out_tmp(2:3);
            clear model_out_tmp
        else
            model_out_above = MLE(errdeg_above,model); %fits without plotting
            model_out_below = MLE(errdeg_below,model); %fits without plotting
        end

        g_aboveN1(i_part,ii) = model_out_above(1);
        g_belowN1(i_part,ii) = model_out_below(1);
        sd_aboveN1(i_part,ii) = model_out_above(2);
        sd_belowN1(i_part,ii) = model_out_below(2);
        
        errdeg_aboveN1{i_part,ii} = errdeg_above;
        errdeg_belowN1{i_part,ii} = errdeg_below;
        
        clear model_out_above model_out_below cutoff above_mdn below_mdn...
            erp_tmp i_elect errdeg_above errdeg_below
    end
    clear ii
end
clear i_part model model2 time_win timelim


% Save variables to existing .mat file
save([saveLocation 'mdn_ERP_split_v4.mat'],'errdeg_aboveN1','errdeg_belowN1','g_aboveN1',...
    'g_belowN1','sd_aboveN1','sd_belowN1','-append')


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% P2 ERP
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

time_win = [200 255]; %P2
timelim = find(times_erp>time_win(1),1)-1:find(times_erp>time_win(2),1)-1;

errdeg_aboveP2 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
errdeg_belowP2 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
g_aboveP2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
g_belowP2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_aboveP2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_belowP2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for selecting specific electrodes
        
        erp_tmp = squeeze(mean(erp_out_all{1,i_part}(i_elect,timelim,:),2)); %average over time
        cutoff = median(erp_tmp); %get median to split trials
        
        above_mdn = find(erp_tmp > cutoff);
        below_mdn = find(erp_tmp <= cutoff);
        
        errdeg_aboveP2{i_part,ii} = resp_errdeg{i_part}(above_mdn);
        errdeg_belowP2{i_part,ii} = resp_errdeg{i_part}(below_mdn);
        
        if i_part == 21 || i_part == 22 %subjects had bias in mean errors
            model_out_tmp = MLE(errdeg_aboveP2{i_part,ii}(:),model2); %fits without plotting
            model_out_aboveP2 = model_out_tmp(2:3);
            clear model_out_tmp
            
            model_out_tmp = MLE(errdeg_belowP2{i_part,ii}(:),model2); %fits without plotting
            model_out_belowP2 = model_out_tmp(2:3);
            clear model_out_tmp
        else
            model_out_aboveP2 = MLE(errdeg_aboveP2{i_part,ii}(:),model); %fits without plotting
            model_out_belowP2 = MLE(errdeg_belowP2{i_part,ii}(:),model); %fits without plotting
        end

        g_aboveP2(i_part,ii) = model_out_aboveP2(1);
        g_belowP2(i_part,ii) = model_out_belowP2(1);
        sd_aboveP2(i_part,ii) = model_out_aboveP2(2);
        sd_belowP2(i_part,ii) = model_out_belowP2(2);
        
        clear model_out_aboveP2 model_out_belowP2 cutoff above_mdn below_mdn...
            erp_tmp i_elect
    end
    clear ii
end
clear i_part model model2 time_win timelim


% Save variables to existing .mat file
save([saveLocation 'mdn_ERP_split_v4.mat'],'errdeg_aboveP2','errdeg_belowP2','g_aboveP2',...
    'g_belowP2','sd_aboveP2','sd_belowP2','-append')


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% N2 ERP
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

time_win = [255 360]; %N2
timelim = find(times_erp>time_win(1),1)-1:find(times_erp>time_win(2),1)-1;

errdeg_aboveN2 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
errdeg_belowN2 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
g_aboveN2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
g_belowN2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_aboveN2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_belowN2 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for selecting specific electrodes
        
        erp_tmp = squeeze(mean(erp_out_all{1,i_part}(i_elect,timelim,:),2)); %average over time
        cutoff = median(erp_tmp); %get median to split trials
        
        above_mdn = find(erp_tmp > cutoff);
        below_mdn = find(erp_tmp <= cutoff);
        
        errdeg_aboveN2{i_part,ii} = resp_errdeg{i_part}(above_mdn);
        errdeg_belowN2{i_part,ii} = resp_errdeg{i_part}(below_mdn);
        
        if i_part == 21 || i_part == 22 %subjects had bias in mean errors
            model_out_tmp = MLE(errdeg_aboveN2{i_part,ii}(:),model2); %fits without plotting
            model_out_aboveN2 = model_out_tmp(2:3);
            clear model_out_tmp
            
            model_out_tmp = MLE(errdeg_belowN2{i_part,ii}(:),model2); %fits without plotting
            model_out_belowN2 = model_out_tmp(2:3);
            clear model_out_tmp
        else
            model_out_aboveN2 = MLE(errdeg_aboveN2{i_part,ii}(:),model); %fits without plotting
            model_out_belowN2 = MLE(errdeg_belowN2{i_part,ii}(:),model); %fits without plotting
        end

        g_aboveN2(i_part,ii) = model_out_aboveN2(1);
        g_belowN2(i_part,ii) = model_out_belowN2(1);
        sd_aboveN2(i_part,ii) = model_out_aboveN2(2);
        sd_belowN2(i_part,ii) = model_out_belowN2(2);
        
        clear model_out_aboveN2 model_out_belowN2 cutoff above_mdn below_mdn...
            erp_tmp i_elect
    end
    clear ii
end
clear i_part model model2 time_win timelim


% Save variables to existing .mat file
save([saveLocation 'mdn_ERP_split_v4.mat'],'errdeg_aboveN2','errdeg_belowN2','g_aboveN2',...
    'g_belowN2','sd_aboveN2','sd_belowN2','-append')



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% P3 ERP
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

time_win = [360 500]; %P3
timelim = find(times_erp>time_win(1),1)-1:find(times_erp>time_win(2),1)-1;

errdeg_aboveP3 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
errdeg_belowP3 = cell(length(exp.participants),length(exp.electrode)); %pre-allocate
g_aboveP3 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
g_belowP3 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_aboveP3 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
sd_belowP3 = NaN(length(exp.participants),length(exp.electrode)); %pre-allocate
model = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); %model with mu
for i_part = 1:length(exp.participants)
    
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for selecting specific electrodes
        
        erp_tmp = squeeze(mean(erp_out_all{1,i_part}(i_elect,timelim,:),2)); %average over time
        cutoff = median(erp_tmp); %get median to split trials
        
        above_mdn = find(erp_tmp > cutoff);
        below_mdn = find(erp_tmp <= cutoff);
        
        errdeg_aboveP3{i_part,ii} = resp_errdeg{i_part}(above_mdn);
        errdeg_belowP3{i_part,ii} = resp_errdeg{i_part}(below_mdn);
        
        if i_part == 21 || i_part == 22 %subjects had bias in mean errors
            model_out_tmp = MLE(errdeg_aboveP3{i_part,ii}(:),model2); %fits without plotting
            model_out_aboveP3 = model_out_tmp(2:3);
            clear model_out_tmp
            
            model_out_tmp = MLE(errdeg_belowP3{i_part,ii}(:),model2); %fits without plotting
            model_out_belowP3 = model_out_tmp(2:3);
            clear model_out_tmp
        else
            model_out_aboveP3 = MLE(errdeg_aboveP3{i_part,ii}(:),model); %fits without plotting
            model_out_belowP3 = MLE(errdeg_belowP3{i_part,ii}(:),model); %fits without plotting
        end

        g_aboveP3(i_part,ii) = model_out_aboveP3(1);
        g_belowP3(i_part,ii) = model_out_belowP3(1);
        sd_aboveP3(i_part,ii) = model_out_aboveP3(2);
        sd_belowP3(i_part,ii) = model_out_belowP3(2);
        
        clear model_out_aboveP3 model_out_belowP3 cutoff above_mdn below_mdn...
            erp_tmp i_elect
    end
    clear ii
end
clear i_part model model2 time_win timelim


% Save variables to existing .mat file
save([saveLocation 'mdn_ERP_split_v4.mat'],'errdeg_aboveP3','errdeg_belowP3','g_aboveP3',...
    'g_belowP3','sd_aboveP3','sd_belowP3','-append')



% /////////////////////////////////////////////////////////////////////////
% Clear workspace
ccc
% /////////////////////////////////////////////////////////////////////////






