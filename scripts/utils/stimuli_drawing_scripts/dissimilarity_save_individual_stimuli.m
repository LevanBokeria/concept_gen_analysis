function dissimilarity_save_individual_stimuli(curr_space,subj,ses,debugMode)

%% Description

% This script will save all 16 of the stimuli as separate PNG files.
% I then used those PNG files to create one giant picture with all 16
% stimuli randomly distributed. This image is intended to be used at the
% beginning of the dissimilarity task, so participants have an overview of
% what types of stimuli to expect.

%% Dissimilarity rating script
% sca;
% clearvars -except ptb_window windowRect; clc; close all;
dbstop if error;

exptStart = tic;
%% Define global variables

% Check you're in the right directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/unified_space/');
end

addpath(genpath(home));

% Define which space you're in, subject ID, and session
if (nargin < 1)
    
    subj = '100';
    ses  = '1';
%     curr_space = 'shade_base_space';
        curr_space = 'beak_tail_space';
    %     curr_space = 'neck_legs_space';
%         curr_space = 'square_circle_space';
%         curr_space = 'pot_leaf_space';
%         curr_space = 'handle_circle_space';
    
    debugMode = 0;
    
end

% Create the appropriate folders
% task = 'dissimilarity';

% Flags

saveImg = 1;

%% Load the appropriate subject and condition specific files

% Make save folders
saveLoc = fullfile(home,'inputs','pre_exposure_imgs',curr_space);

if ~exist(saveLoc)
    mkdir(saveLoc);
end

[pairs, pair_indices, coords, nSesPairs] = genDissim_full;


%% Start PTB

% Setup PTB with some default values
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2);

if debugMode
    PsychDebugWindowConfiguration(0,0.5);
end

ptb_window = 10;
windowRect = [0 0 1280 1024];

if numel(Screen('Screens')) > 1
    whichScreen = 1;
else
    whichScreen = max(Screen('Screens'));
end

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
%         [ptb_window, windowRect] = PsychImaging('OpenWindow', whichScreen, white, [], 32, 2,...
%             [], [],  kPsychNeed32BPCFloat);
            [ptb_window, windowRect] = PsychImaging('OpenWindow', whichScreen, white);

end

if ~strcmp(curr_space,'handle_circle_space')
    Screen('BlendFunction', ptb_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end

hz = Screen('NominalFrameRate',ptb_window);

if ~debugMode
    HideCursor;
end

%--- Get the KbQueue up and running.
KbQueueCreate();
KbQueueStart();

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', whichScreen);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% make screen white
Screen(ptb_window,'FillRect',white);
Screen('Flip',ptb_window);

% set keys
spaceKey  = KbName('space');
escapeKey = KbName('ESCAPE'); %exit game

% Show or hide the mouse
ShowCursor('Arrow');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SPACE SPECIFIC PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define instructions based on space
switch curr_space
    
    case 'shade_base_space'
        
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_shade_base;
        
        % Stimulus parameters
        [base_width, shade_height, ...
            switch_trunk_offset,...
            trunk_width,trunk_height,...
            switch_width,switch_height,...
            knob_diameter,overlap] ...
            = parameters_shade_base(SF_all,SF_body);
        
        % LOAD ALL PARTS OF THE lamp
        [base,base_idx,shade,shade_idx,Bar_Neck_idx, Bar_Legs_idx] ...
            = load_shade_base_parts(home,ptb_window);
        
        stimName = 'lamp';
        bar_x_name = 'base';
        bar_y_name = 'shade';
        
        % Location of the lamps on the screen
        nOptions = 2;
        
        target_stim_offset_x = 250 * SF_all; % distance between the testing and target lamps in x dimension
        target_stim_offset_y = 150; % distance between the testing and target lamps in y dimension
        
        x_offset = (screenXpixels - (base_width * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
        y_offset = (screenYpixels - (target_stim_offset_y + base_width)*2) / 2;
        
        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions);
        %    op_y_loc = repmat(yCenter,1,nOptions) + 328 * SF_all; %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;
        
        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] - 100;
        
        
        if nOptions == 1
            op_x_loc = xCenter;
            op_y_loc = yCenter;
        end
      
    case 'beak_tail_space'
        
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

        stimName   = 'bird';
        bar_x_name = 'beak';
        bar_y_name = 'tail';
        
        % Location of the birds on the screen
        nOptions = 2;
        
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
        
    case 'flower_space'
        
       % Enable alpha blending so transparency works
        Screen('BlendFunction', ptb_window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        % Get the flower spaces
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_flower;
        
        % Parameters of the flower
        [Stem_width, StemTip_height, StemTip_width, StemTip_overlap, StemTip_offset, ...
            flower_height, flower_width, ...
            Stem_body_offset, overlap] ...
            = parameters_flower(SF_all,SF_body);
        
        % Load flower parts and make them into a texture
        [Stem, StemTip, Heart, Body, Bar_Stem, Bar_Heart, ...
            Stem_idx, StemTip_idx, Heart_idx, Body_idx, ...
            Bar_x_idx, Bar_y_idx] ...
            = load_flower_parts(home,ptb_window);        

        stimName   = 'flower';
        bar_x_name = 'stem';
        bar_y_name = 'heart';
        
        % Determine centroid display location of 3 options.
        nOptions = 2;
        
        target_flower_offset_x = 100; % distance between the testing and target flowers in x dimension
        target_flower_offset_y = 0; % distance between the testing and target flowers in y dimension
        
        x_offset = (screenXpixels - (flower_width * (nOptions-1) + target_flower_offset_x)) / 2;
        y_offset = (screenYpixels - (target_flower_offset_y + flower_height)*2) / 2;
        
        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions) + -20;
        op_y_loc = repmat(yCenter,1,nOptions); %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;
        
        if nOptions == 1
            op_x_loc = xCenter;
            op_y_loc = yCenter;
        end
        
    case 'neck_legs_space'

        % Enable alpha blending so transparency works
        Screen('BlendFunction', ptb_window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        % Get the bird spaces
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_neck_legs;

        % Parameters of the bird
        [bird_height,bird_width, ...
            neck_width, legs_width, ...
            feet_width, feet_height, ...
            head_height, head_width, ...
            neck_body_offset, ...
            head_neck_offset, head_neck_overlap, ...
            legs_body_offset, legs_feet_offset, ...
            overlap] ...
            = parameters_neck_legs(SF_all,SF_body);


        % Load bird parts and make them into a texture
        [Head, Neck, Body, Legs, Feet, ...
            Head_idx, Neck_idx, Legs_idx, Body_idx, Feet_idx, ...
            Bar_x_idx, Bar_y_idx] ...
            = load_neck_legs_parts(home,ptb_window);
               
        stimName   = 'bird';
        bar_x_name = 'neck';
        bar_y_name = 'legs';              
        
        % Location of the birds on the screen
        nOptions = 2;

        target_stim_offset_x = 50; % distance between the testing and target flowers in x dimension
        target_stim_offset_y = 300 * SF_all; % distance between the testing and target flowers in y dimension

        x_offset = (screenXpixels - (bird_width * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
        y_offset = (screenYpixels - (target_stim_offset_y + bird_height)*2) / 2;

        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions) - 40;
        op_y_loc = repmat(yCenter + 60 ,1,nOptions); %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;

        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2];


        if nOptions == 1
            op_x_loc = xCenter;
            op_y_loc = yCenter;
        end
        
        % Border line parameters
        border_length = 150;
        lineWidthPix  = 2;
        
        cussion = 50;
        cussion2 = 40;
        
        xCoords =  [-border_length + cussion2, border_length + cussion2, ...
            -border_length + cussion2, border_length + cussion2];
        yCoords =  [-(bird_height/2 + max(scaled_space_x) + cussion), -(bird_height/2 + max(scaled_space_x) + cussion), ...
            bird_height/2 + max(scaled_space_x) + cussion, bird_height/2 + max(scaled_space_x) + cussion];
        allCoords = [xCoords;yCoords];
        
    case 'square_circle_space'

        
        [Bar_x_idx, Bar_y_idx] ...
            = load_square_circle_parts(ptb_window);
        
        % Gen forage parameters
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_square_circle;
        
        % Gen object parameters
        [square_diameter, coeff, baseRect] = parameters_square_circle;
       
        % Location of the squircles on the screen
        nOptions = 2;
        
        target_stim_offset_x = 150 * SF_all; % distance between the testing and target flowers in x dimension
        target_stim_offset_y = 150; % distance between the testing and target flowers in y dimension
        
        x_offset = (screenXpixels - (square_diameter * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
        y_offset = (screenYpixels - (target_stim_offset_y + square_diameter)*2) / 2;
        
        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions);
        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] - 80;
        %     op_y_loc = repmat(yCenter,1,nOptions) - 150; %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;
        
        if nOptions == 1
            op_x_loc = xCenter;
            op_y_loc = yCenter - 150;
        end
        
        circle_x_loc = op_x_loc;
        circle_y_loc = op_y_loc + square_diameter*coeff/2 + 30*SF_all;
        
        % Names
        stimName = 'squircle';
        bar_x_name = 'square';
        bar_y_name = 'circle';       
        
    case 'pot_leaf_space'
        
        [flower,leaf,flower_idx,leaf_idx, ...
            Bar_x_idx, Bar_y_idx] ...
            = load_pot_leaf_parts(home,ptb_window);
        
        % Gen forage parameters
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_pot_leaf;
        
        % Gen object parameters
        [pot_bottom,pot_height,...
            trunk_width,trunk_height,trunk_rectColor,...
            stem_width,stem_height,stem_rectColor,...
            flower_width,flower_height,...
            leaf_width,overlap] ...
            = parameters_pot_leaf(SF_all,SF_body);
        
        % Location of the birds on the screen
        nOptions = 2;

        target_stim_offset_x = 550 * SF_all; % distance between the testing and target flowers in x dimension
        target_stim_offset_y = 150; % distance between the testing and target flowers in y dimension

        x_offset = (screenXpixels - (pot_bottom * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
        y_offset = (screenYpixels - (target_stim_offset_y + pot_bottom)*2) / 2;

        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions);
        op_y_loc = repmat(yCenter,1,nOptions) + 328 * SF_all; %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;

        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] + 100;


        if nOptions == 1
            op_x_loc = xCenter;
            op_y_loc = yCenter + 380;
        end

        % Names
        stimName   = 'plant';
        bar_x_name = 'pot size';
        bar_y_name = 'leaf size';
        
    case 'handle_circle_space'
        % Space parameters
        [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
            stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
            max_navig, lower_edge_idx, upper_edge_idx] ...
            = genForageParams_handle_circle;
        
        
        % Stimulus parameters
        [handle_height,handle_width,handleBaseRect,rod_height,rod_width,rodBaseRect,...
            handleFrameColor,rodColor,circleColor,penWidthPixels,...
            sigma,contrast,aspectRatio,phase,orientation,gabortex] = parameters_handle_circle(ptb_window);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Location of the squircles on the screen
        nOptions = 1;
        
        target_stim_offset_x = 0 * SF_all; % distance between the testing and target flowers in x dimension
        target_stim_offset_y = 250; % distance between the testing and target flowers in y dimension
        
        x_offset = (screenXpixels - (handle_width/2 * (nOptions-1) + target_stim_offset_x * (nOptions - 1))) / 2;
        y_offset = (screenYpixels - (target_stim_offset_y + handle_height)*2) / 2;
        
        op_x_loc = linspace(x_offset,screenXpixels - x_offset,nOptions) - 100;
        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2];
        %     op_y_loc = repmat(yCenter,1,nOptions) - 150; %   linspace(y_offset,screenYpixels - y_offset,nOptions) + -90;
        
        if nOptions == 1
            op_x_loc = xCenter - 130;
            op_y_loc = yCenter;
        end
        
        circle_x_loc = op_x_loc + handle_width/2 + rod_width;
        circle_y_loc = op_y_loc;
        
end % switch

% on_screen_instructions_line = ['How dissimilar are these ' stimName 's overall?'];

% Create parameters for the scale
scale_max_val = 10;

switch curr_space
    case 'shade_base_space'
        scaleY = repmat(op_y_loc(end) + scaled_space_y(250) + 170,1,scale_max_val);
    case 'beak_tail_space'
        scaleY = repmat(op_y_loc(end) + bird_height/2 + 100,1,scale_max_val);            
    case 'flower_space'
        scaleY = repmat(op_y_loc(end) - flower_height/2 - 100,1,scale_max_val);
    case 'neck_legs_space'
        scaleY = repmat(op_y_loc(end) - bird_height - 500,1,scale_max_val);
    case 'square_circle_space'
        scaleY = repmat(op_y_loc(end) - square_diameter,1,scale_max_val);  
    case 'pot_leaf_space'
        scaleY = repmat(op_y_loc(end) + pot_bottom + 100,1,scale_max_val);               
end
scaleX = linspace(250,screenXpixels -250,scale_max_val);
button_size = 20;

% Define location of the red dot
baseRect_mouse = [0 0 30 30];

% For Ovals we set a miximum diameter up to which it is perfect for
maxDiameter = max(baseRect_mouse);

%% Run main trials

for iStim = 1:size(coords,1)
    
    % Display the two flowers
    xPos1 = coords(iStim,1);
    yPos1 = coords(iStim,2);
    
    % Temp
    %                     xPos1 = max(coords(:));
    %                     yPos1 = max(coords(:));
    
    xPos1_val = scaled_space_x(round(xPos1));
    yPos1_val = scaled_space_y(round(yPos1));
    
    % Display the stimuli
    switch curr_space
        
        case 'shade_base_space'
            
            % Draw auto
            [trunk_rect_centered,base_rect_centered,...
                shade_rect_centered,switch_rect_centered,...
                knob_rect_centered] = draw_shade_base(...
                SF_all,SF_body,...
                xPos1_val,yPos1_val,op_x_loc(1),op_y_loc(1));
            
            % Draw everything
            Screen('FillRect', ptb_window, [0 0 0], trunk_rect_centered);
            Screen('DrawTextures',ptb_window, [shade_idx,base_idx],[], ...
                [shade_rect_centered;base_rect_centered]');
            
            % Draw the rect to the screen
            Screen('FillRect', ptb_window, [0 0 0], switch_rect_centered);
            Screen('FillOval', ptb_window, [0 0 0], knob_rect_centered);
        
        case 'beak_tail_space'
            
            [rBeak, rBody, rTail, rBeakTip, rTailTip] ...
                = draw_beak_tail(SF_all,SF_body,...
                xPos1_val,yPos1_val,...
                op_x_loc(1),op_y_loc(1));
            
            Screen('DrawTextures',ptb_window,[Body_idx, Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx], [], [rBody; rBeak; rTail; rBeakTip; rTailTip]');
            
        case 'flower_space'
            
            [rStem, rStemTip, rBody, rHeart] ...
                = drawFlower(SF_all,SF_body,...
                xPos1_val,yPos1_val, ...
                op_x_loc(1),op_y_loc(1));
            
            Screen('DrawTextures', ptb_window, [Stem_idx, StemTip_idx Body_idx, Heart_idx], [], [rStem; rStemTip; rBody; rHeart]');
            
        case 'neck_legs_space'
            
            [rBody, rNeck, rHead, rLegs, rFeet] ...
                = draw_neck_legs(SF_all,SF_body,...
                xPos1_val,yPos1_val,...
                op_x_loc(1),op_y_loc(1));
            Screen('DrawTextures', ptb_window, [Head_idx, Neck_idx Body_idx Legs_idx Feet_idx], [], [rHead; rNeck; rBody; rLegs; rFeet]');
                        
        case 'square_circle_space'
            
            % Draw the test object
            [rectColor,centeredRect, ovalColor, centeredOval] = ...
                draw_square_circle(xPos1_val,yPos1_val,op_x_loc(1), op_y_loc(1),...
                circle_x_loc(1),circle_y_loc(1));
            
            Screen('FillRect', ptb_window, rectColor, centeredRect);
            Screen('FillOval', ptb_window, ovalColor, centeredOval);
  
        case 'pot_leaf_space'
            
            % Test
            draw_pot_leaf(ptb_window,SF_all,SF_body,...
                flower_idx,leaf_idx,...
                xPos1_val,yPos1_val,op_x_loc(1), op_y_loc(1));
            
        case 'handle_circle_space'
            
            [propertiesMat,centeredHandle,centeredRod,...
                rectColor,baseOval,centeredOval] = ...
                draw_handle_circle(ptb_window,xPos1_val,yPos1_val,...
                op_x_loc(1),op_y_loc(1),...
                circle_x_loc(1),circle_y_loc(1));
            
    end % switch
    
    pause(0.1);
    %% Save a screenshot of the trials

    switch curr_space
        
        case 'shade_base_space'
            
            % GetImage call. Alter the rect argument to change the location of the screen shot
%             rect_to_save = [op_x_loc(1) - square_diameter/2 - 5, ...
%                             op_y_loc(1) - square_diameter*coeff/2 - 5, ...
%                             op_x_loc(1) + square_diameter/2 + 5, ...
%                             circle_y_loc(1) + scaled_space_y(250)/2 + 5];
%                         
                        
            rect_to_save = [op_x_loc(1) - scaled_space_x(250)/2 - 5, ...
                op_y_loc(1) - trunk_height/2 - shade_height - 5, ...
                op_x_loc(1) + scaled_space_x(250)/2 + 5, ...
                op_y_loc(1) + trunk_height/2 + scaled_space_y(250) + 5];

            
        case 'pot_leaf_space'
            
            % GetImage call. Alter the rect argument to change the location of the screen shot
            rect_to_save = [op_x_loc(1) - scaled_space_x(250)/2 + overlap, ...
                op_y_loc(1) - pot_height/2 - trunk_height - flower_height - overlap, ...
                op_x_loc(1) + scaled_space_x(250)/2 + overlap, ...
                op_y_loc(1) + pot_height/2];
            
        case 'neck_legs_space'
            
            % GetImage call. Alter the rect argument to change the location of the screen shot
            rect_to_save = [op_x_loc(1) - bird_width/2 + 5, ...
                op_y_loc(1) - bird_height/2 - head_height - scaled_space_x(250) - overlap - 3, ...
                op_x_loc(1) + bird_width/2 + overlap + 75, ...
                op_y_loc(1) + feet_height + bird_height/2 + scaled_space_x(250) - 8];
                        
        case 'beak_tail_space'

            % GetImage call. Alter the rect argument to change the location of the screen shot
            rect_to_save = [op_x_loc(1) - bird_width/2 - beak_body_offset_x - scaled_space_x(250) + overlap - beakTip_width - 5, ...
                op_y_loc(1) - bird_height/2 - 5, ...
                op_x_loc(1) + bird_width/2 - overlap + scaled_space_y(250) + tailTip_width + 5, ...
                op_y_loc(1) + bird_height/2 + 5];
                        
        case 'square_circle_space'
            % GetImage call. Alter the rect argument to change the location of the screen shot
            rect_to_save = [op_x_loc(1) - square_diameter/2 - 5, ...
                op_y_loc(1) - square_diameter*coeff/2 - 5, ...
                op_x_loc(1) + square_diameter/2 + 5, ...
                circle_y_loc(1) + scaled_space_y(250)/2 + 5];
            
        case 'handle_circle_space'
            
            rect_to_save = [op_x_loc(1) - handle_width/2 - 5, ...
                circle_y_loc(1) - scaled_space_y(250)/2 - 5,...
                circle_x_loc(1) + scaled_space_y(250)/2 + 5,...
                circle_y_loc(1) + scaled_space_y(250)/2 + 5];
            
    end

    imageArray = Screen('GetImage', ptb_window, rect_to_save, 'backBuffer');
    
    Screen('Flip',ptb_window);

    
    if saveImg
        % imwrite is a Matlab function, not a PTB-3 function
        imwrite(imageArray, fullfile(saveLoc,['stim_' int2str(iStim) '.png']));
    end
    pause(0.1);
end % iStims
% KbWait 
sca;

%% Extra functions

    function drawExtraStuff(current_pairs,xPos1,yPos1,xPos2,yPos2,...
            ptb_window,pairDisplayed,nPairs_inner,black,on_screen_instructions_line,...
            idx_i_real_dist,idx_inner)
        
        %     if ~exist('idx_inner','var')
        %         idx_inner = 0;
        %     end
        %
        %     Screen('TextSize',ptb_window,15);
        %
        %     % Draw all text
        %     DrawFormattedText(ptb_window,['Trial ' int2str(pairDisplayed) ' out of ' int2str(nPairs_inner)],20,50);
        %
        %     % Draw correct response
        % %     DrawFormattedText(ptb_window,['Correct vs registered ' ...
        % %         int2str(idx_i_real_dist) ' - ' int2str(idx_inner)],xCenter - 100,yCenter + 250);
        %
        %     if ~practice
        %         DrawFormattedText(ptb_window,['Session ' ses ' out of ' int2str(nSessions_total)],20,80);
        %     end
        % %         DrawFormattedText(ptb_window,int2str(pairDisplayed), 20, 20, black);
        %     DrawFormattedText(ptb_window,on_screen_instructions_line, 'center', 30, black);
        %
        % %     % Draw nTrials per PA
        % %     formatSpec = ['Current pairs %d %d \n\n xPos1 %d \n\n yPos1 %d \n\n xPos2 %d \n\n yPos2 %d'];
        % %     DrawFormattedText(ptb_window,...
        % %         sprintf(formatSpec,[current_pairs(1),current_pairs(2),...
        % %         xPos1,yPos1,xPos2,yPos2]),...
        % %         50, screenYpixels - 150);
        
        
    end

end