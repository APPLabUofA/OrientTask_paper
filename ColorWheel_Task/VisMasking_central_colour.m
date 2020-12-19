%% VisEntrainment_quest.m 
% Dr. Kyle Mathewson - University of Alberta
% Entrainment with repetitive annulus stimulation. 



clear all
close all
% warning off MATLAB:DeprecatedLogicalAPI
% Screen('Preference', 'SkipSyncTests', 1); %for test on an LCD
Priority(2);

% Input participant's number, name and date of testing
% Output file will be named after the inputs

Info.number = 0; %input('Participant Number:','s');
Info.date = datestr(now,30); % 'dd-mmm-yyyy HH:MM:SS' 

% Create text file to write data in
Filename = [Info.number '--' Info.date '_data.mat'];

%% Variables to adjust
%target
vis_thresh = 15; %target RGB points darker than grey background
targ_length = 1; %in refresh cycles; each refresh is 1000msec / 120 Hz = 8.333 msec
target_size = 18;

%mask
mask_thresh = 127; %points darker than background
SOA = 1; %refresh target onset to mask onset (50 ms optimal/8.3333 msec = 5 cycles 
ISI = SOA - targ_length; %target OFFSET to mask onset 
mask_length = 3; %refresh cycles of mask
entrainer_size = 35;

%fixations
fixation_length = 48; %400 ms
preblank_length = 24; %200 ms

%trial numbers    
n_trials = 10; %the +1 accounts for 1 catch trial for every type of target trial
n_blocks = 1; %how many blocks (1st block is a practice session 


%% Open the main window and get dimensions
white=[255,255,255];  %WhiteIndex(window);
black=[0,0,0];   %BlackIndex(window);

grey= [128 128 128]; %background colour
mask_grey = grey-mask_thresh; %mask colour
targ_grey = grey-vis_thresh; %stim_grey - vis_thresh; %target colour

%load the window up
screenNumber = max(Screen('Screens')); % Get the maximum screen number i.e. get an external screen if avaliable
[window,rect]=Screen(screenNumber ,'OpenWindow',grey(1));
HideCursor; %Comment this out for testing    

v_res = rect(4);
h_res = rect(3);
v_center = v_res/2;
h_center = h_res/2;
fixation = [h_center-10 v_center-10];
trigger_size = [0 0 1 1]; %use [0 0 1 1] for eeg, 100 100 for photodiode

% Get presentation timing information
refresh = Screen('GetFlipInterval', window); % Get flip refresh rate
slack = refresh/2; % Divide by 2 to get slack

%% Set up the offScreen windows and stimuli

left_arrow = [h_center-10 v_center-30; h_center+10 v_center-40; h_center+10 v_center-20];   %create the attention cues
right_arrow = [h_center+10 v_center-30; h_center-10 v_center-40; h_center-10 v_center-20];

%setup the blank screen
blank=Screen(window,'OpenoffScreenWindow',grey);
    
%setup the fixation screen
fixate=Screen(window,'OpenoffScreenWindow',grey);
    Screen(fixate, 'DrawLines', [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,

%Setup the mask
mask=Screen(window,'OpenoffScreenWindow',grey);
    Screen(mask, 'FillOval',[mask_grey(1) grey(1); mask_grey(2) grey(2); mask_grey(3) grey(3)] , [h_center-entrainer_size h_center-target_size; v_center-entrainer_size v_center-target_size; h_center+entrainer_size h_center+target_size; v_center+entrainer_size v_center+target_size] );
    
    %make the target screen for the actual task
target=Screen(window,'OpenoffScreenWindow',grey);
    Screen(target, 'FillOval',[1 0 0], [h_center-target_size v_center-target_size h_center+target_size v_center+target_size] ); %Target
    disp(h_res);
    disp(v_res);
%     disp(v_center-target_size);
%     disp(h_center+target_size);
%     disp(v_center+target_size);

    
WaitSecs(2);

%% Instructions
Screen('DrawLines',window, [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,
Screen('DrawText',window,'Please keep your eyes fixed on the central cross the entire time.',fixation(1)-500,fixation(2)-100,0);  %Display instructions
Screen('Flip', window,[],0); %flip it to the screen
 KbWait; %wait for subject to press button

%wait a bit before it starts
Screen('CopyWindow',blank ,window,rect,rect);
Screen(window, 'Flip')
WaitSecs(2);

%% Loop for blocks
for i_block = 1:n_blocks;
    
   
    %% Loop for trials
    for i_trial = 1:n_trials 
        
   
        %% put up fixation
        Screen('CopyWindow',fixate ,window,rect,rect);
        tfixate_onset = Screen(window, 'Flip');
        
        %% put up a blank
        Screen('CopyWindow',blank ,window,rect,rect);
        tblank_onset = Screen(window, 'Flip', tfixate_onset + fixation_length*refresh - slack);
        fixation_time(i_block,i_trial) = tblank_onset - tfixate_onset; 
        

        %% present the Target
        Screen('CopyWindow',target ,window,rect,rect);
        ttarget_onset = Screen(window, 'Flip', tblank_onset + preblank_length*refresh - slack);
        blank_time(i_block,i_trial) = ttarget_onset - tblank_onset;
        Screen('FillRect',window,[1 0 0],trigger_size);
        

         %% blank Inter stimulus interval
        Screen('CopyWindow',blank ,window,rect,rect);
        tISI_onset = Screen(window, 'Flip', ttarget_onset + targ_length*refresh - slack);
        target_time(i_block,i_trial) = tISI_onset - ttarget_onset;

        %% present the mask
        Screen('CopyWindow',mask ,window,rect,rect);
        tmask_onset = Screen(window, 'Flip', tISI_onset + ISI*refresh - slack);
        ISI_time(i_block,i_trial) = tmask_onset - tISI_onset;
        
        %% Response period
        Screen('CopyWindow',blank ,window,rect,rect);
        Screen(window, 'Flip', tmask_onset + mask_length*refresh - slack);
       
        t1 = GetSecs;
        keyIsDown = 0;
        while  ~keyIsDown 
            [keyIsDown, secs, keyCode] = KbCheck;
        end
        response = find(keyCode>0);   %1 is 49, 5 is 53, left arrow is 37, right is 39
        
        %keep a log of the subject answers
        if response == 37
            subject_answer(i_block,i_trial) = 1; %detected
        elseif response == 39
            subject_answer(i_block,i_trial) = 0; %undetected
        else
            subject_answer(i_block,i_trial) = 99;
        end
        subject_rt(i_block,i_trial) = secs-t1; %compute response time and log
        WaitSecs(1); %wait one second before next trial
        
    end
    

    %% Wait for the subject to move onto the next block, or end the experiment.   
    if i_block == 1;    
        Screen('DrawLines',window, [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,
        Screen('DrawText',window,'Let the experimenter know if you have any questions. Press any key to begin the next practise block of the experiment.',fixation(1)-800,fixation(2)-110,0);  
        Screen('Flip', window,[],0); %
        WaitSecs(1);
        KbWait 
    elseif i_block == 2;    
        Screen('DrawLines',window, [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,
        Screen('DrawText',window,'Let the experimenter know if you have any questions. Press any key to begin the first block of the experiment.',fixation(1)-800,fixation(2)-110,0);  
        Screen('Flip', window,[],0); %
        WaitSecs(1);
        KbWait 
    elseif i_block == n_blocks
        Screen('DrawLines',window, [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,
        Screen('DrawText',window,['You have now completed all ' num2str(n_blocks-2) ' blocks. Press any key to end the experiment.'],fixation(1)-600,fixation(2)-110,0);  
        Screen('Flip', window,[],0); %
        WaitSecs(1);
        KbWait 
    else
        Screen('DrawLines',window, [-7 7 0 0; 0 0 -7 7], 1, 0, [h_center,v_center],0);  %Print the fixation,
        if i_block-1 == 1
            suffix = '';
        else
            suffix = 's';
        end
        Screen('DrawText',window,['You have now completed ' num2str(i_block-1) ' block' suffix ' out of ' num2str(n_blocks-1) '. Press any key to start the next block.'],fixation(1)-600,fixation(2)-110,0);  
        Screen('Flip', window,[],0); %
        WaitSecs(1);
        KbWait 
        
    end
    %wait a bit before it goes on
    Screen('CopyWindow',blank ,window,rect,rect);
    Screen(window, 'Flip')
    WaitSecs(2);

end




%% Save the data and close the window, exit
Screen('Close', window);
ShowCursor;
save(Filename,'subject_answer','subject_lags','subject_present','subject_rt','Info', 't', 'target_time', 'entr_time', 'entrblank_time', 'lag_time', 'fixation_time');

%compute average RT and detection and plot
for i_lag = 1:4; 
    mean_per_lag(i_lag) =  sum(subject_answer((subject_lags==i_lag)&subject_present==1)==1)/length(subject_answer((subject_lags==i_lag)&subject_present==1)); 
    RT_per_lag(i_lag) = mean(subject_rt((subject_lags==i_lag)&subject_present==1));
end
figure; subplot(2,1,1); plot(1-mean_per_lag); ylabel('MissProportion');
        subplot(2,1,2); plot(RT_per_lag*1000); ylabel('RT (msec)');
 

