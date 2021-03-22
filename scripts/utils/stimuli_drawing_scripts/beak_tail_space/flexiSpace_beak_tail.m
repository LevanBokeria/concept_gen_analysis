%% Explore the display of the bird.

% Script allows to play with different beak tail sizes, body sizes, and
% global scaling. To see what fits best on the screens. 
% sca
clear; clc; close all;

sesStart = tic;

% Check you're in the right directory 
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/unified_space/');
end

addpath(genpath(home));

%% Setup PTB
% Setup PTB with some default values
PsychDefaultSetup(2);   

Screen('Preference','SkipSyncTests',2);

debugMode = 1;

if debugMode    
    PsychDebugWindowConfiguration(0,0.8);
end

ptb_window = 10;
windowRect = [0 0 1280 1024];

if numel(Screen('Screens')) > 1
    whichScreen = 2;
else
    whichScreen = max(Screen('Screens'));
end
% 
try
    info = Screen('GetWindowInfo',ptb_window);
    screenOpen = 1;
catch
    screenOpen = 0;
end

% set colors
black = BlackIndex(whichScreen); % pixel value for black
white = WhiteIndex(whichScreen);

% screen size
if screenOpen == 0
    [ptb_window, windowRect] = PsychImaging('OpenWindow', whichScreen, white);
end

Screen('BlendFunction', ptb_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', whichScreen);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

Screen('Flip',ptb_window);

%% Task parameters
automatic = 1;

%% If automatic
if automatic

    % Get the bird spaces
    [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
        stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
        max_navig, lower_edge_idx, upper_edge_idx] ...
        = genForageParams_beak_tail;


    % Parameters of the bird
    [beak_height, beakTip_width, beakTip_height, ...
              tail_height, tailTip_width, ...
              bird_height, bird_width, ...
              beak_body_offset_y, beak_body_offset_x, ...
              tail_body_offset, beakTip_offset, tailTip_offset, ...
              overlap] ...
              = parameters_beak_tail(SF_all,SF_body);

    % Load bird parts and make them into a texture
    [Beak, BeakTip, Tail, TailTip, Body, Bar_Beak, Bar_Tail, ...
        Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx, Body_idx, ...
        Bar_x_idx, Bar_y_idx] ...
        = load_beak_tail_parts(home,ptb_window);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Location of the birds on the screen
    nOptions = 1;

    target_stim_offset_x = -200 * SF_all; % distance between the testing and target flowers in x dimension
    target_stim_offset_y = 350; % distance between the testing and target flowers in y dimension

    x_offset = (screenXpixels - (bird_width * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
    y_offset = (screenYpixels - (target_stim_offset_y + bird_width)*2) / 2;

    op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions);
    op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] + 0;
    
    
    if nOptions == 1
        op_x_loc = xCenter;
        op_y_loc = yCenter;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    step_x = 6;
    
    mult1x = 10 - step_x;
    mult2x = 9  - step_x;
    
    step_y = step_x;
    
    mult1y = 10 - step_y;
    mult2y = 9  - step_y;
    
    xPos = 25*mult1x;
    yPos = 25*mult1y;
    
    xPos_val = scaled_space_x(xPos);
    yPos_val = scaled_space_y(yPos);
    
    % xPos_val = max(scaled_space_x);
    % yPos_val = max(scaled_space_y);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % Draw all options
    for iOption = 1:nOptions
          
        [rBeak, rBody, rTail, rBeakTip, rTailTip] ...
        = draw_beak_tail(SF_all,SF_body,...
                   xPos_val,yPos_val,...
                   op_x_loc(iOption),op_y_loc(iOption));

        Screen('DrawTextures',ptb_window,[Body_idx, Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx], [], [rBody; rBeak; rTail; rBeakTip; rTailTip]');

    end % iOption

%     % Remaining script
%     % Dials
%     bar_x_name = 'pot';
%     bar_y_name = 'leaves';
%     
%     draw_controller_decoupled(...
%         ptb_window,...
%         bar_x_name,bar_y_name,...
%         Bar_x_idx,Bar_y_idx,...
%         x_bar_pos,y_bar_pos,...
%         bar_limits,...
%         black)
%     
    Screen('Flip',ptb_window);      
    
    
%% If manual
else
    
    %% Forage parameters
    sf_mult = -6;

    SF_all  = 1 + sf_mult/10; % scaling factor for everything together, including the dimensions.
    SF_body = 3; % scaling factor for the size of the body, excluding the range of dimensions and minimum sizes of the dimensions.
    SF_dim  = 1; % scaling factor for dimensions
    
    %set parameters
    rotationAngle     = 0; %the initial rotation angle in radians
    xPos              = 0; %coordinates of the bird on the top
    yPos              = 0;
    delta_translation = 1; %increment in translation. This used to be multiplied by SF_dim, but now that we're standardizing the space and instead scaling the
    delta_rotation    = 0.04; %increment in rotation
    dist_criterion    = 25; % distance around the reward in the standard space.
    
    % Linear scale on both dimensions 
    % Generate the standardized space
    nDiscrete = 300;

    % x dimension (beak size)

    % we want the 50th value to be the minimum value. and 300th to be the
    % maximum
    xPos_min = 25;
    xPos_max = 800;

    stand_space_x = linspace(xPos_min,xPos_max,nDiscrete);
    scaled_space_x = stand_space_x * SF_all; % .* linspace(1,1.6,nDiscrete);

%     stand_space_x  = linspace(xPos_min,xPos_max, nDiscrete + xPos_min);
%     scaled_space_x = stand_space_x;

    % y dimension (tail size)
    yPos_min = xPos_min;
    yPos_max = xPos_max;

    stand_space_y = linspace(yPos_min,yPos_max, nDiscrete);

    scaled_space_y = stand_space_y * SF_all; % .* linspace(1,1.6,nDiscrete);

    % Define the navigable space (i.e. boundaries).
    % These will be indices to look up the corresponding number in the
    % scaled_space vector
    max_navig = 300;

    lower_edge_idx = 50;
    upper_edge_idx = 300;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    beak_height    = 24 * SF_all * SF_body; % 50 pixesl
    beakTip_width  = 100 * SF_all * SF_body;
    beakTip_height = 42 * SF_all * SF_body;
    
    tail_height   = beak_height;
    tailTip_width = beakTip_width;
    tailTip_height = beakTip_height;
    
    bird_height = 214 * SF_all * SF_body;
    bird_width  = 150 * SF_all * SF_body;
    
    beak_body_offset_y = 16 * SF_all * SF_body;
    beak_body_offset_x = -5 * SF_all * SF_body;
    tail_body_offset = 128 * SF_all * SF_body;
    beakTip_offset   = 1 * SF_all * SF_body;
    tailTip_offset   = -1 * SF_all * SF_body;
    
    overlap = 1; % overlap between bird parts, to avoid the 1 pixel white line
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Location of the birds on the screen
    nOptions = 1;

    target_stim_offset_x = -200 * SF_all; % distance between the testing and target flowers in x dimension
    target_stim_offset_y = 350; % distance between the testing and target flowers in y dimension

    x_offset = (screenXpixels - (bird_width * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
    y_offset = (screenYpixels - (target_stim_offset_y + bird_width)*2) / 2;

    op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions);
    op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] + 0;
    
    
    if nOptions == 1
        op_x_loc = xCenter;
        op_y_loc = yCenter + 200;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    step_x = 8;
    
    mult1x = 10 - step_x;
    mult2x = 9  - step_x;
    
    step_y = step_x;
    
    mult1y = 10 - step_y;
    mult2y = 9  - step_y;
    
    xPos = 25*mult1x;
    yPos = 25*mult1y;
    
    xPos_val = scaled_space_x(xPos);
    yPos_val = scaled_space_y(yPos);
    
    % xPos_val = max(scaled_space_x);
    % yPos_val = max(scaled_space_y);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load bird parts and make them into a texture
    [Beak, BeakTip, Tail, TailTip, Body, Bar_Beak, Bar_Tail, ...
        Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx, Body_idx, ...
        Bar_x_idx, Bar_y_idx] ...
        = load_beak_tail_parts(home,ptb_window);
    
    % Draw all options
    for iOption = 1:nOptions

        % Define the boxes
        rBody = [0, ...
                 0, ...
                 bird_width, ...
                 bird_height];

        rBody = CenterRectOnPoint(rBody, op_x_loc(iOption),op_y_loc(iOption));

        rBeak = [rBody(1) + overlap - xPos_val - beak_body_offset_x, ...
                 rBody(2) + beak_body_offset_y, ...
                 rBody(1) + overlap - beak_body_offset_x, ...
                 rBody(2) + beak_body_offset_y + beak_height];

        rBeakTip = [rBeak(1) - beakTip_width + overlap, ...
                    rBeak(2) + beakTip_offset, ...
                    rBeak(1) + overlap, ...
                    rBeak(2) + beakTip_offset + beakTip_height];

        rTail = [rBody(3), ...
                 rBody(2) + tail_body_offset, ...
                 rBody(3) + yPos_val, ...
                 rBody(2) + tail_body_offset + tail_height];

        rTailTip = [rTail(3) - overlap, ...
                    rTail(4) - beakTip_height + tailTip_offset, ...
                    rTail(3) - overlap + tailTip_width, ...
                    rTail(4) + tailTip_offset];           

        Screen('DrawTextures',ptb_window,[Body_idx, Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx], [], [rBody; rBeak; rTail; rBeakTip; rTailTip]');

    end % iOption

    Screen('Flip', ptb_window);    
    
    
end % if automatic



