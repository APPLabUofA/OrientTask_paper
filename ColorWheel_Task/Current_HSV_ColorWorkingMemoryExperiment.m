% COLORWORKINGMEMORYEXPERIMENT Runs a color working memory task
% a la Zhang & Luck (2008). The task requires memory for the color of
% briefly presented squares. Participants then report the color of a single
% probed square using a continuous report task.
%
%   ColorWorkingMemoryExperiment();
%
%	Preferences can be found down at the bottom, beginning on line 197.
%



function ColorWorkingMemoryExperiment()

colour_test;

  try
    prepareEnvironment;
    window = openWindow();
    prefs = getPreferences();
    
     % Get presentation timing information
     refresh = Screen('GetFlipInterval', window.onScreen); % Get flip refresh rate
     slack = refresh/2; % Divide by 2 to get slack

    % Put up instructions and wait for keypress.
    instruct(window);
    returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize)

    jitter1 = randi([0,5])*0.1;
    
   WaitSecs(0.5 + jitter1);

    % Get rects for each item.
    rects = cell(1, max(prefs.setSizes));
    for i = 1:max(prefs.setSizes)
      rects{i} = circularArrayRects([0, 0, prefs.squareSize, prefs.squareSize], ...
        i, prefs.radius, window.centerX, window.centerY)';
    end

    colorWheelLocations = colorwheelLocations(window,prefs);

    for trialIndex = 1:length(prefs.fullFactorialDesign)
        

      % Determine how many items there are on this trial and the duration.
      nItems = prefs.setSizes(prefs.fullFactorialDesign(prefs.order(trialIndex), 1));
      retentionInterval = prefs.retentionIntervals; %(prefs.fullFactorialDesign(prefs.order(trialIndex), 2));
      SOA = prefs.SOA(prefs.fullFactorialDesign(prefs.order(trialIndex), 2));
      ISI = SOA - prefs.targ_length; %target OFFSET to mask onset

      % Pick an item to test.
      itemToTest(trialIndex) = RandSample(1:nItems);

      % Pick the colors for this trial.
      colorsInDegrees{trialIndex} = ceil(rand(1, nItems)*360);
      
      % Draw fixation.
%         drawFixation(window, window.centerX, window.centerY, prefs.fixationSize);
        
         %Interval
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
       tblank_onset = Screen('Flip', window.onScreen, prefs.fixation_length*refresh - slack);
        

 if SOA == 0
      % Draw the stimulus.
       colorsToDisplay = prefs.colorwheel(colorsInDegrees{trialIndex}, :)';
       Screen('FillOval', window.onScreen, [0 0 0], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth]);
       Screen('FillOval', window.onScreen, colorsToDisplay, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(1),prefs.trigger_size);
       
      % Post the stimulus 
       ttarget_onset = Screen('Flip', window.onScreen, tblank_onset + prefs.preblank_length*refresh - slack);
           
       %No Interval
     
       %Mask
       Screen('FillOval', window.onScreen, [0 0 0], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth]);
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);

       t_maskonset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);
      
        %Interval
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
       Screen('Flip', window.onScreen, t_maskonset + (prefs.mask_length - prefs.targ_length)*refresh - slack);
       
       WaitSecs(retentionInterval);
       end
       
       if SOA > 0
           % Draw the stimulus.
       colorsToDisplay = prefs.colorwheel(colorsInDegrees{trialIndex}, :)';
       Screen('FillOval', window.onScreen, colorsToDisplay, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(SOA),prefs.trigger_size);
       
      % Post the stimulus 
       ttarget_onset = Screen('Flip', window.onScreen, tblank_onset + prefs.preblank_length*refresh - slack);
           
         %Interval
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
       tISI_onset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);


       %Mask
       Screen('FillOval', window.onScreen, [0 0 0], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth]);
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);

       t_maskonset = Screen('Flip', window.onScreen, tISI_onset + ISI*refresh - slack);   
       
       %Interval
       Screen('FillOval', window.onScreen, window.gray, rects{nItems});
       Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
       Screen('Flip', window.onScreen, t_maskonset + prefs.mask_length*refresh - slack);
       
       WaitSecs(retentionInterval);
       end
      % Remove stimulus, return to blank, wait for retention interval to pass.
%        returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize);
       
%        

      % Choose a circle to test, then display the response screen.
      data.presentedColor(trialIndex) = deg2rad(colorsInDegrees{trialIndex}(itemToTest(trialIndex)));
      colorsOfTest = repmat([120 120 120], nItems, 1);
      colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
%       drawFixation(window, window.centerX, window.centerY, prefs.fixationSize);
      Screen('FillRect', window.onScreen, colorsOfTest', rects{nItems});

      drawColorWheel(window, prefs);

      % Wait for click.
      SetMouse(window.centerX, window.centerY);
      ShowCursor('Arrow');

      % If mouse button is already down, wait for release.
      GetMouse(window.onScreen);
      buttons = 0;
      while any(buttons)
        [x, y, buttons] = GetMouse(window.onScreen);
      end

      everMovedFromCenter = false;
      while ~any(buttons)

        drawColorWheel(window, prefs);

        [x,y,buttons] = GetMouse(window.onScreen);
        [minDistance, minDistanceIndex] = min(sqrt((colorWheelLocations(1, :) - x).^2 + (colorWheelLocations(2, :) - y).^2));

        if(minDistance < 250)
          everMovedFromCenter = true;
        end

        if(everMovedFromCenter)
          colorsOfTest(itemToTest(trialIndex), :) = prefs.colorwheel(minDistanceIndex,:);
        else
          colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
        end

%         drawFixation(window, window.centerX, window.centerY, prefs.fixationSize);
        Screen('FillOval', window.onScreen, colorsOfTest', rects{nItems});
        
        drawColorWheel(window, prefs);
        Screen('Flip', window.onScreen);
        
      end
      data.reportedColor(trialIndex) = deg2rad(minDistanceIndex);
      while any(buttons) % wait for release
        [x,y,buttons] = GetMouse(window.onScreen);
      end

      HideCursor
      
      if trialIndex == 50
            rest(window);
        end
      
        if trialIndex == 100
            rest(window);
        end
        
        if trialIndex == 150
            rest(window);
        end
        
        if trialIndex == 200
            rest(window);
        end
        
        if trialIndex == 250
            rest(window);
        end
        
        if trialIndex == 300
            rest(window);
        end
        
        if trialIndex == 350
            rest(window);
        end

      % Return to fixation.
      returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize);
      
      jitter2 = randi([0,5])*0.1;
   
      WaitSecs(0.5 + jitter2);
    end

    % Preliminary analysis of results.
    data.error = (180/pi) .* (angle(exp(1i*data.reportedColor)./exp(1i*data.presentedColor)));
    data.setSize = prefs.setSizes(prefs.fullFactorialDesign(prefs.order, 1));
    data.retentionInterval = prefs.retentionIntervals; %(prefs.fullFactorialDesign(prefs.order,2));
    data.SOA = prefs.SOA(prefs.fullFactorialDesign(prefs.order,2));
    
    save data.mat data prefs
    postpareEnvironment;

  catch
    postpareEnvironment;
    psychrethrow(psychlasterror);
    
  end % end try/catch
end % end whole colorworkingmemoryscript

function prepareEnvironment

  clear all;
  HideCursor;

  commandwindow; % Select the command window to avoid typing in open scripts

  % Seed the random number generator.
  RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));

%   ListenChar(2); % Don't print to MATLAB command window
%   Screen('Preference', 'SkipSyncTests', 1);
%   
       %% Set up parallel port
    %initialize the inpoutx64 low-level I/O driver
    config_io;
    %optional step: verify that the inpoutx64 driver was successfully installed
    global cogent;
    if( cogent.io.status ~= 0 )
       error('inp/outp installation failed');
    end
    %write a value to the default LPT1 printer output port (at 0x378)
    address_eeg = hex2dec('B010');

    outp(address_eeg,0);  %set pins to zero  
end

function postpareEnvironment
  ShowCursor;
%   ListenChar(0);
  Screen('CloseAll');
end

function instruct(window)
  prefs = getPreferences(); 
  Screen('TextSize', window.onScreen, window.fontsize);
  Screen('DrawText', window.onScreen, 'Remember the colors. Click to begin.', 100, 100, 255);
  Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
  Screen('Flip', window.onScreen);
  GetClicks(window.onScreen);
end

function rest(window)
  prefs = getPreferences();
  Screen('TextSize', window.onScreen, window.fontsize);
  Screen('DrawText', window.onScreen, 'Feel free to take a break at this time.  Click the mouse when you are ready to continue.', 100, 100, 255);
  Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
  Screen('Flip', window.onScreen);
  GetClicks(window.onScreen);
end

function drawFixation(window, fixationX, fixationY, fixationSize)
prefs = getPreferences();
Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(92), prefs.trigger_size);
end

function offsets = circularArrayOffsets(n, centerX, centerY, radius, rotation)
  degreeStep = 360/n;
  offsets = [sind(0:degreeStep:(360-degreeStep) + rotation)'.* radius, ...
             cosd(0:degreeStep:(360-degreeStep) + rotation)'.* radius];
end

function rects = circularArrayRects(rect, nItems, radius, centerX, centerY)
  coor = circularArrayOffsets(nItems, centerX, centerY, radius, 0) + repmat([centerX centerY], nItems, 1);
  rects = [coor(:, 1)-rect(3)/2 , coor(:, 2)-rect(3)/2, coor(:, 1)+rect(3)/2, coor(:, 2)+rect(3)/2];
end

function returnToFixation(window, fixationX, fixationY, fixationSize)
  prefs = getPreferences();
  Screen('FillRect', window.onScreen, window.gray);
  Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
  Screen('FillRect',window.onScreen, Vpixx2Vamp(92), prefs.trigger_size);
  Screen('Flip', window.onScreen);
end

function window = openWindow()
  window.screenNumber = max(Screen('Screens'));
  window.onScreen = Screen('OpenWindow', window.screenNumber, [125 125 125]);
  Screen('BlendFunction', window.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
  [window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
  window.screenRect  = [0, 0, window.screenX, window.screenY]; % screen rect
  window.centerX = window.screenX * 0.5; % center of screen in X direction
  window.centerY = window.screenY * 0.5; % center of screen in Y direction
  window.centerXL = floor(mean([0, window.centerX])); % center of left half of screen in X direction
  window.centerXR = floor(mean([window.centerX, window.screenX])); % center of right half of screen in X direction

  % Basic drawing and screen variables.
  window.black    = BlackIndex(window.onScreen);
  window.white    = WhiteIndex(window.onScreen);
  window.gray     = mean([window.black window.white])
  window.fontsize = 24;
  window.bcolor   = window.gray;
end

function drawColorWheel(window, prefs)
prefs = getPreferences();
  colorWheelLocations = [cosd(1:360).*prefs.colorWheelRadius + window.centerX; ...
    sind(1:360).*prefs.colorWheelRadius + window.centerY];
  colorWheelSizes = 20;
  Screen('DrawDots', window.onScreen, colorWheelLocations, colorWheelSizes, prefs.colorwheel', [], 1);
  Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
end

function L = colorwheelLocations(window,prefs)
  L = [cosd(1:360).*prefs.colorWheelRadius + window.centerX; ...
       sind(1:360).*prefs.colorWheelRadius + window.centerY];
end

function prefs = getPreferences
    prefs.nTrialsPerCondition = 135;
    prefs.setSizes = [1];
    prefs.retentionIntervals = [0.50];
    prefs.squareSize = 20; % size of each stimulus object, in pixels
    prefs.radius = 0;
    prefs.fixationSize = 2;
    prefs.fixation_length = 120;

    %target
    prefs.targ_length = 1; %in refresh cycles; each refresh is 1000msec / 120 Hz = 8.333 msec

    %mask
    prefs.SOA = [0, 4, 10]; %refresh target onset to mask onset (50 ms optimal/8.3333 msec = 5 cycles
    prefs.maskwidth = 15;
    prefs.mask_length = 6; %refresh cycles of mask
    % Variables

    prefs.trigger_size = [0 0 1 1];
    
    %fixations
    prefs.preblank_length = 24; %200 ms
    
    % Colorwheel details.
    prefs.colorWheelRadius = 350;
    prefs.colorwheel = load('newcolorwheel360.mat', 'HSV_wheel');
    prefs.colorwheel = prefs.colorwheel.HSV_wheel;

    % Randomize trial order of full factorial design order.
    prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
        length(prefs.SOA), ...
        prefs.nTrialsPerCondition]);

    prefs.order = Shuffle(1:length(prefs.fullFactorialDesign));
    prefs.nTrials = length(prefs.fullFactorialDesign);
end