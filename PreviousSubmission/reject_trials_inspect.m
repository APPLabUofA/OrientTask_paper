function EEG = reject_trials_inspect(EEG,part_name,exp,rejtrial,i_set)

% Additional rejection of trials for participants after a visual review of
% the trials using pop_eegplot(EEG) then keeping a record of the trial
% numbers I selected for rejection. The trial number goes into the 
% EEG.reject.rejthresh(#)=1; function where # is replaced by the trial
% number. 
% 
% Need to set breakpoint at function in Preprocessing_OrientTask.m to be 
% able to make changes to this code without exiting the function

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Select participant and processing settings to remove list of trials

if (strcmpi(part_name,'subj3') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(16)=1;
    EEG.reject.rejthresh(100)=1;
    EEG.reject.rejthresh(167)=1;
    EEG.reject.rejthresh(169)=1;
    EEG.reject.rejthresh(173)=1;
    EEG.reject.rejthresh(234)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1); %creates a list of trials rejected
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); %actually rejects the trials
    
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------    
    
elseif (strcmpi(part_name,'subj5') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(52)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj5') && strcmpi(exp.settings,'filt_byTargets_v4'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(172)=1;
    EEG.reject.rejthresh(185)=1;
    EEG.reject.rejthresh(194)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
    
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj6') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(19)=1;
    EEG.reject.rejthresh(58)=1;
    EEG.reject.rejthresh(61)=1;
    EEG.reject.rejthresh(68)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj6') && strcmpi(exp.settings,'filt_byTargets_v4'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(70)=1;
    EEG.reject.rejthresh(83)=1;
    EEG.reject.rejthresh(113:115)=1;
    EEG.reject.rejthresh(125)=1;
    EEG.reject.rejthresh(136)=1;
    EEG.reject.rejthresh(163)=1;
    EEG.reject.rejthresh(191)=1;
    EEG.reject.rejthresh(195)=1;
    EEG.reject.rejthresh(201)=1;
    EEG.reject.rejthresh(207:208)=1;
    EEG.reject.rejthresh(218)=1;
    EEG.reject.rejthresh(220:221)=1;
    EEG.reject.rejthresh(238)=1;
    EEG.reject.rejthresh(241)=1;
    EEG.reject.rejthresh(244)=1;
    EEG.reject.rejthresh(261:263)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
    
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj7') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(39)=1;
    EEG.reject.rejthresh(50:51)=1;
    EEG.reject.rejthresh(64:72)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj9') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(152)=1;
    EEG.reject.rejthresh(155)=1;
    EEG.reject.rejthresh(234)=1;
    EEG.reject.rejthresh(238)=1;
    EEG.reject.rejthresh(261)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj10') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(48)=1;
    EEG.reject.rejthresh(67)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj10') && strcmpi(exp.settings,'filt_byTargets_v4'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(30)=1;
    EEG.reject.rejthresh(79)=1;
    EEG.reject.rejthresh(191:193)=1;
    EEG.reject.rejthresh(250)=1;
    EEG.reject.rejthresh(265)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj11') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(48:49)=1;
    EEG.reject.rejthresh(115)=1;
    EEG.reject.rejthresh(121)=1;
    EEG.reject.rejthresh(182)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj12') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(47)=1;
    EEG.reject.rejthresh(68)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj12') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(48:49)=1;
    EEG.reject.rejthresh(64)=1;
    EEG.reject.rejthresh(96)=1;
    EEG.reject.rejthresh(108)=1;
    EEG.reject.rejthresh(118)=1;
    EEG.reject.rejthresh(246)=1;
    EEG.reject.rejthresh(250:252)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj13') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(33)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj13') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(10:11)=1;
    EEG.reject.rejthresh(30)=1;
    EEG.reject.rejthresh(151:152)=1;
    EEG.reject.rejthresh(192)=1;
    EEG.reject.rejthresh(230:231)=1;
    EEG.reject.rejthresh(234)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj14') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(51)=1;
    EEG.reject.rejthresh(55:56)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj14') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(91)=1;
    EEG.reject.rejthresh(174)=1;
    EEG.reject.rejthresh(189:190)=1;
    EEG.reject.rejthresh(201)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj15') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(57)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj15') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(58:59)=1;
    EEG.reject.rejthresh(81:83)=1;
    EEG.reject.rejthresh(104)=1;
    EEG.reject.rejthresh(149)=1;
    EEG.reject.rejthresh(151)=1;
    EEG.reject.rejthresh(153:154)=1;
    EEG.reject.rejthresh(160:161)=1;
    EEG.reject.rejthresh(197)=1;
    EEG.reject.rejthresh(200:203)=1;
    EEG.reject.rejthresh(205:206)=1;
    EEG.reject.rejthresh(228)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj16') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(43)=1;
    EEG.reject.rejthresh(53)=1;
    EEG.reject.rejthresh(120)=1;
    EEG.reject.rejthresh(230)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj17') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(55)=1;
    EEG.reject.rejthresh(63)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj17') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(2)=1;
    EEG.reject.rejthresh(47)=1;
    EEG.reject.rejthresh(66)=1;
    EEG.reject.rejthresh(101)=1;
    EEG.reject.rejthresh(152)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj18') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(70)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj18') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(101)=1;
    EEG.reject.rejthresh(257)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(31)=1;
    EEG.reject.rejthresh(82:85)=1;
    EEG.reject.rejthresh(93)=1;
    EEG.reject.rejthresh(104:108)=1;
    EEG.reject.rejthresh(115:116)=1;
    EEG.reject.rejthresh(139)=1;
    EEG.reject.rejthresh(146)=1;
    EEG.reject.rejthresh(154:158)=1;
    EEG.reject.rejthresh(181)=1;
    EEG.reject.rejthresh(186:188)=1;
    EEG.reject.rejthresh(191)=1;
    EEG.reject.rejthresh(195:199)=1;
    EEG.reject.rejthresh(220:221)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(14)=1;
    EEG.reject.rejthresh(20)=1;
    EEG.reject.rejthresh(25:26)=1;
    EEG.reject.rejthresh(28)=1;
    EEG.reject.rejthresh(41:42)=1;
    EEG.reject.rejthresh(46)=1;
    EEG.reject.rejthresh(58)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj21') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(47)=1;
    EEG.reject.rejthresh(69)=1;
    EEG.reject.rejthresh(71)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj21') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(30:31)=1;
    EEG.reject.rejthresh(116)=1;
    EEG.reject.rejthresh(118)=1;
    EEG.reject.rejthresh(145:147)=1;
    EEG.reject.rejthresh(175)=1;
    EEG.reject.rejthresh(179:181)=1;
    EEG.reject.rejthresh(192)=1;
    EEG.reject.rejthresh(206:207)=1;
    EEG.reject.rejthresh(212:213)=1;
    EEG.reject.rejthresh(231)=1;
    EEG.reject.rejthresh(239)=1;
    EEG.reject.rejthresh(249)=1;
    EEG.reject.rejthresh(251)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj22') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(4)=1;
    EEG.reject.rejthresh(40)=1;
    EEG.reject.rejthresh(95)=1;
    EEG.reject.rejthresh(99)=1;
    EEG.reject.rejthresh(104)=1;
    EEG.reject.rejthresh(153)=1;
    EEG.reject.rejthresh(245)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav'))%only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(9)=1;
    EEG.reject.rejthresh(18)=1;
    EEG.reject.rejthresh(24)=1;
    EEG.reject.rejthresh(33)=1;
    EEG.reject.rejthresh(40)=1;
    EEG.reject.rejthresh(60)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    %subject moved throughout experiment
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(10)=1;
    EEG.reject.rejthresh(23)=1;
    EEG.reject.rejthresh(26)=1;
    EEG.reject.rejthresh(39)=1;
    EEG.reject.rejthresh(42)=1;
    EEG.reject.rejthresh(51)=1;
    EEG.reject.rejthresh(57)=1;
    EEG.reject.rejthresh(72:73)=1;
    EEG.reject.rejthresh(87)=1;
    EEG.reject.rejthresh(107)=1;
    EEG.reject.rejthresh(110)=1;
    EEG.reject.rejthresh(117)=1;
    EEG.reject.rejthresh(132:133)=1;
    EEG.reject.rejthresh(141)=1;
    EEG.reject.rejthresh(149)=1;
    EEG.reject.rejthresh(156)=1;
    EEG.reject.rejthresh(158)=1;
    EEG.reject.rejthresh(167)=1;
    EEG.reject.rejthresh(170:172)=1;
    EEG.reject.rejthresh(174)=1;
    EEG.reject.rejthresh(178)=1;
    EEG.reject.rejthresh(180:182)=1;
    EEG.reject.rejthresh(188:189)=1;
    EEG.reject.rejthresh(194)=1;
    EEG.reject.rejthresh(201)=1;
    EEG.reject.rejthresh(214)=1;
    EEG.reject.rejthresh(220)=1;
    EEG.reject.rejthresh(228)=1;
    EEG.reject.rejthresh(231)=1;
    EEG.reject.rejthresh(241)=1;
    EEG.reject.rejthresh(246)=1;
    EEG.reject.rejthresh(252:253)=1;
    EEG.reject.rejthresh(258)=1;
    EEG.reject.rejthresh(263:264)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj24') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(69)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj24') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(119)=1;
    EEG.reject.rejthresh(192)=1;
    EEG.reject.rejthresh(228)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj25') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    %subject wore a lot of makeup
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(13)=1;
    EEG.reject.rejthresh(17)=1;
    EEG.reject.rejthresh(19)=1;
    EEG.reject.rejthresh(70)=1;
    EEG.reject.rejthresh(98)=1;
    EEG.reject.rejthresh(101)=1;
    EEG.reject.rejthresh(208)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj26') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(14:15)=1;
    EEG.reject.rejthresh(53)=1;
    EEG.reject.rejthresh(115)=1;
    EEG.reject.rejthresh(240)=1;
    EEG.reject.rejthresh(269)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj27') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(2)=1;
    EEG.reject.rejthresh(7)=1;
    EEG.reject.rejthresh(75)=1;
    EEG.reject.rejthresh(112)=1;
    EEG.reject.rejthresh(118)=1;
    EEG.reject.rejthresh(135)=1;
    EEG.reject.rejthresh(140:141)=1;
    EEG.reject.rejthresh(143)=1;
    EEG.reject.rejthresh(170:171)=1;
    EEG.reject.rejthresh(186:188)=1;
    EEG.reject.rejthresh(190)=1;
    EEG.reject.rejthresh(229)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj27') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(43)=1;
    EEG.reject.rejthresh(48:49)=1;
    EEG.reject.rejthresh(51)=1;
    EEG.reject.rejthresh(61)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);     

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj29') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(1)=1;
    EEG.reject.rejthresh(52)=1;
    EEG.reject.rejthresh(65)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj30') && strcmpi(exp.settings,'filt_byTargets_v4')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(113:114)=1;
    EEG.reject.rejthresh(211:212)=1;
    EEG.reject.rejthresh(232)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------

elseif (strcmpi(part_name,'subj30') && strcmpi(exp.settings,'filt_byCatchTargets_v2_wav')) %only when using this specific settings 
    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
    EEG.reject.rejthresh(65)=1;
    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
end %elseif loop end

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%% Save rejected trials
EEG.rejtrial = rejtrial;

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
