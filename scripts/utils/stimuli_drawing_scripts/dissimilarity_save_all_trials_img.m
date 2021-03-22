function dissimilarity_save_all_trials_img(curr_space,subj,ses,debugMode)

%% Description

% This function will generate PNG images of the pairs of stimuli used in
% the dissimilarity/similarity experiments.

% It can generate both the practice trial and real trial images.

% There are flags for whether to include the correct response in the image
% that is being generated.


%% Dissimilarity rating script
% sca;
clearvars -except ptb_window windowRect; clc; close all;
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
%     curr_space = 'square_circle_space';
%     curr_space = 'pot_leaf_space';
%     curr_space = 'handle_circle_space';
        
    
    debugMode = 0;
    
end    

% Flags
print_corr_resp     = 0; % prints what the correct response is on the screen
run_practice_trials = 0; % whether or not to print practice trials
run_actual_trials   = 1; % whether or not to print the real trials
save_data           = 1; % flag for debuggning. if 0 then stuff doesn't get saved.

% Create the appropriate folders 
task = 'dissimilarity';

nSessions_total = 4;
nReps = nSessions_total/2;

% Create dissimilarity matrices 
dm = cell(1,nReps);

for iCell = 1:nReps
    dm{iCell} = zeros(16,16);
end

%% Load the appropriate subject and condition specific files

% Make save folders
saveLoc = fullfile(home,'inputs','stimuli_screenshots',curr_space);

if ~exist(saveLoc) %#ok<*EXIST>
    mkdir(saveLoc);
end
              
if ~exist(fullfile(saveLoc,'parameters'))
    mkdir(fullfile(saveLoc,'parameters'));
    mkdir(fullfile(saveLoc,'seeds'));    
end

% Save the seed.
seed = rng('shuffle'); % to reproduce
if save_data
    save(fullfile(saveLoc,'seeds',['seed_' task ses '.mat']),'seed');
end

% Load the task table
% dissimilarity_table = load(fullfile(home,'inputs','pregenerated_files',task,[task '_tables.mat']));
% dissimilarity_table = dissimilarity_table.dissimilarity_table;
% 
% dissimilarity_table = dissimilarity_table.congruent{1,phase};
% 
% pairs     = dissimilarity_table.pairs;
% coords    = dissimilarity_table.coords;
% nSesPairs = dissimilarity_table.nSesPairs;

% if ~strcmp(ses,'1')
%     % So subject has run session 1, then load the parameters so we can
%     % continue from where it stopped
% 
%     load(fullfile(saveLoc,'parameters','dissimilarity1.mat'),'pairs','pair_indices','coords','nSesPairs');
%     
% else
%     % Generate the trial sequence, randomized for every different subject
% %     [pairs, pair_indices, coords, nSesPairs] = genDissim_full;
%     
%     % Temporary to repeat what happened in pot_leaf_space
%     load('C:\Users\lb08\Documents\GitHub\con_learn\unified_space\outputs\100\pot_leaf_space\dissimilarity\parameters\dissimilarity1.mat',...
%         'pairs','pair_indices','coords','nSesPairs');
% 
%     
% end
% Load the example trials
% example_trial_pairs = load(fullfile(home,'inputs','pregenerated_files',task,'example_trials.mat'));
% example_trial_pairs = example_trial_pairs.example_trials;

% Create pairs and coords and nSesPair, for saving all pairs including ID
% pairs.

pairs = load(fullfile(home,'\inputs\pre_generated_files\dissimilarity_pairs.mat'))
pairs = pairs.pairs;

pairs1 = pairs(1:120,:);
pairs2 = pairs(121:240,:);

% Sanity check
pairs1_sorted = unique(pairs1,'rows');
pairs2_flip = fliplr(pairs2);
pairs2_sorted = unique(pairs2_flip,'rows');
assert(isequal(pairs1_sorted,pairs2_sorted));

pairs1 = [pairs1; repmat([1:16]',1,2)];
pairs2 = [pairs2; repmat([1:16]',1,2)];

pairs = [pairs1; pairs2];

pair_indices = [1:size(pairs,1)]';

coords = [100,100;100,150;100,200;100,250;150,100;150,150;150,200;150,250;200,100;200,150;200,200;200,250;250,100;250,150;250,200;250,250];

nSesPairs = size(pairs,1) / nSessions_total;

%% Creat matrices and files for initial exposure 
exemplarDur = 1; %#ok<*NASGU>
ITI         = 0;   

% Trial duration
trial_dur = 6;
keyIsDown = 0;
%% Start PTB

% Setup PTB with some default values
PsychDefaultSetup(2);      
Screen('Preference', 'SkipSyncTests', 2);

if debugMode
    PsychDebugWindowConfiguration(0,0.5);
end

ptb_window = 10;
windowRect = [0 0 1280 720];

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
            = genForageParams_square_circle; %#ok<*ASGLU>
        
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

        op_y_loc = [yCenter - target_stim_offset_y/2, yCenter + target_stim_offset_y/2] + 150;


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
        nOptions = 2;
        
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
        
        % Names
        stimName   = 'rackets';
        bar_x_name = 'spatial frequency';
        bar_y_name = 'circle size';        
        
        
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
        scaleY = repmat(op_y_loc(end) + square_diameter + scaled_space_x(250),1,scale_max_val);  
    case 'pot_leaf_space'
        scaleY = repmat(op_y_loc(end) + pot_bottom + 100,1,scale_max_val);  
    case 'handle_circle_space'
        scaleY = repmat(op_y_loc(end) + scaled_space_y(250) + 100,1,scale_max_val);  
end
scaleX = linspace(250,screenXpixels -250,scale_max_val);
button_size = 20;        

% Define location of the red dot
baseRect_mouse = [0 0 30 30];

% For Ovals we set a miximum diameter up to which it is perfect for
maxDiameter = max(baseRect_mouse);    

%% Run example trials

if strcmp(ses,'1') && run_practice_trials
    
    practice = 1;
    
    iRep = 1;
    
    % Trial duration
    trial_dur_practice = 1;
    
    missed_trials = 0;
    
    % Edit save location
    saveLoc_practice = fullfile(saveLoc,'example_imgs');
    
    if ~exist(fullfile(saveLoc_practice,'parameters'))
        mkdir(fullfile(saveLoc_practice,'parameters'));
    end
    
    %% Define trials
    
    % Give specific example trials.
    pairs_current = example_trial_pairs;
    
    for iExTrial = 1:size(pairs_current,1)
    
        [matched,match_idx] = ismember(pairs_current(iExTrial,:),pairs,'rows');
        
        if ~matched
            [matched,match_idx] = ismember(fliplr(pairs_current(iExTrial,:)),pairs,'rows');
        end
        
        pair_indices_current(iExTrial) = pair_indices(match_idx); %#ok<*AGROW>
    
    end
    
    nPairs_current = size(pairs_current,1);
    
    %% Run the trials
    
    % Save parameters first
    curr_date_time = datetime;
    if save_data
        save(fullfile(saveLoc_practice,'parameters',[task ses '.mat']));
    end
    
    [responseTable_meta] = run_dissimilarity(saveLoc_practice,ses,pairs_current,...
        pair_indices_current,nPairs_current,missed_trials,trial_dur_practice,practice);

    % Draw text
    Screen('TextSize',ptb_window,12);

    DrawFormattedText(ptb_window,['This is the end of the practice trials!\n\n', ...
        'Remember, now you will have 6 seconds to respond on each trial!',...
        '\n\n You will do total of ' int2str(nSessions_total) ' sessions, during which you will see ' ...
        '120 pairs of stimuli, each repeated 3 times.',...
        '\n\n So total number of trials will be ' int2str(120*3) ...
        '\n\n Press any key to start real trials!'],...
        'center','center');
    
    Screen('Flip',ptb_window);
    KbWait([], 3);
    
    KbReleaseWait();
    
end

%% Run main trials

if run_actual_trials

    practice = 0; %#ok<*UNRCH>

    % Do a session loop

    for iSes = str2double(ses):nSessions_total

        % Define the repetition
        if ismember(iSes,[1,2])
            iRep = 1;
        elseif ismember(iSes,[3,4])
            iRep = 2;
        elseif ismember(iSes,[5,6])
            iRep = 3;
        end

        missed_trials  = 0;

        sesStart = tic;

        responseTable_meta = table;

        % Pairs for the current session
        pairs_current = pairs((iSes-1)*nSesPairs+1:(iSes)*nSesPairs,:);

        pair_indices_current = pair_indices((iSes-1)*nSesPairs+1:(iSes)*nSesPairs);

        nPairs_current = size(pairs_current,1);

        
        %% Run the trials
        
        % Save all the parameters
        curr_date_time = datetime;
        if save_data
            save(fullfile(saveLoc,'parameters',[task ses '.mat']));
        end
        
        [responseTable] = run_dissimilarity(saveLoc,ses,pairs_current,...
            pair_indices_current,nPairs_current,missed_trials,trial_dur,practice);
        
        %% Check if need to repeat
        
        responseTable_meta = [responseTable_meta; responseTable];
        % Adjust the trial column
        responseTable_meta.trial = (1:height(responseTable_meta))';

        %% Present break
        Screen('TextSize',ptb_window,12);

        if iSes == nSessions_total
            % Draw text
            DrawFormattedText(ptb_window,['This is the end of the experiment!.\n\n' 'Thank you!\n\n' 'Press any key to exit'],'center','center');
            Screen('Flip',ptb_window);
            KbWait([], 3);

            KbReleaseWait();
            sca;

        elseif iSes == 2 || iSes == 4

            second_rep_text_2 = ['Please take a 30 second break.' ,...
                '\n\n\n\n For the next 2 sessions, you will be tested once more on the same', ...
                ' 120 pairs of ' stimName 's.',...
                '\n\n\n The ratings you give to the upcoming stimuli do not necessarily have to be consistent ', ...
                '\n\n with the ratings you gave to the same stimuli in the previous sessions.',...
                '\n\n If you start to have a different feeling of what number on the scale corresponds with what degree of perceived dissimilarity,',...
                '\n\n feel free to adjust your responses accordingly. You do not have to give the same responses',...
                '\n\n as you gave before. Always give the response that you feel is correct in the moment.',...
                '\n\n\n For example, if for a certain level of perceived dissimilarity you previously used 10',...
                '\n\n but now you think that 8 is a more appropriate response, you should respond with 8.',...
                '\n\n\n\n Press any button once ready to continue.'];

            DrawFormattedText(ptb_window,second_rep_text_2,'center','center');

%             drawNumScale(ptb_window,scaleX,scaleY,button_size,black);

            Screen('Flip',ptb_window);

            KbReleaseWait;

        else

            break_text = ['This is the end of session ' ses '!', ...
                '\n\n Press any button to start session ' int2str(str2double(ses) + 1) '.'];

            DrawFormattedText(ptb_window,break_text,'center','center');
            Screen('Flip',ptb_window);

    %         KbWait([],3);
            KbReleaseWait;

        end

        % Update ses
        ses = int2str(iSes + 1);


    end % for iSes

end % if save_actual_trials

KbReleaseWait;
sca;

%% Function to run trials
function [responseTable_temp] = run_dissimilarity(saveLoc_inner,ses,...
        pairs_current_inner,pair_indices_current_inner,...
        nPairs_current_inner,missed_trials_inner,trial_dur_inner,practice) %#ok<*INUSL>

        responseTable_temp  = table;

        if practice
            saveFolder = 'example_imgs';
        else
            saveFolder = 'pair_imgs';
        end

        if print_corr_resp
            saveFolder = [saveFolder '_correct_responses'];
        end
        
        if ~practice
            
            % Save the images in the appropriate sub-order folder
            saveFolder = fullfile(saveFolder,['order_' int2str(iRep)]);
            
            % Create that folder if it doesn't exist
            if save_data
                if ~exist( fullfile(saveLoc,saveFolder))
                    mkdir((fullfile(saveLoc,saveFolder)));
                end
            end
        end
        
        for iPair = 1:nPairs_current_inner

            % Set mouse position to the center
            [~, my, ~] = GetMouse(ptb_window);
            SetMouse(xCenter, my, ptb_window);

            % Record trial start time
            trialStartTime = GetSecs;

            % Reset previous trial choice 
            idx = [];

            %% Draw the graph
%                 close all
%                 figure
%                 D = zeros(16,16);
                real_distances = squareform(pdist(coords));
                unique_dist = unique(real_distances);
                
%                 index_mat = [4,8,12,16; 3,7,11,15; 
                
                stim_1 = pairs_current(iPair,1);
                stim_2 = pairs_current(iPair,2);
                
%                 D(stim_1,stim_2) = 1;
%                 D(stim_2,stim_1) = 1;
                
                i_real_dist = real_distances(stim_1,stim_2);
                idx_i_real_dist = find(unique_dist == i_real_dist);
                
                % For similarity scale, reverse it
                reversed_scale = 10:-1:1;
                idx_i_real_dist = reversed_scale(idx_i_real_dist);
                
%                 G = graph(D);
%                 G.Edges.Weight = idx_i_real_dist;
%                 
%                 h = plot(G,'EdgeLabel',G.Edges.Weight);
%                 h.XData = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
%                 h.YData = [repmat([1 2 3 4],1,4)];
%                 
%                 title(int2str(idx_i_real_dist));
%                 
%                 set(gcf,'Position',[  -747   672   560   420]);
%                 pause(0.2);

            %% do the timer and present the stimuli
            while (GetSecs - trialStartTime) < trial_dur_inner


                % Display the two flowers
                xPos1 = coords(pairs_current_inner(iPair,1),1);
                yPos1 = coords(pairs_current_inner(iPair,1),2);

                % Temp
%                     xPos1 = max(coords(:));
%                     yPos1 = max(coords(:));

                xPos1_val = scaled_space_x(round(xPos1));
                yPos1_val = scaled_space_y(round(yPos1));

                xPos2 = coords(pairs_current_inner(iPair,2),1);
                yPos2 = coords(pairs_current_inner(iPair,2),2);

                % Temp
%                     xPos2 = max(coords(:));
%                     yPos2 = max(coords(:));                    

                xPos2_val = scaled_space_x(round(xPos2));
                yPos2_val = scaled_space_y(round(yPos2));

                % Say which pairs are shown
%                 drawExtraStuff(pairs_current_inner(iPair,:),xPos1,yPos1,xPos2,yPos2)
                
                % Draw the scale and text
%                 drawNumScale(ptb_window,scaleX,scaleY,button_size,black);
%                 drawExtraStuff(pairs_current_inner(iPair,:),xPos1,yPos1,xPos2,yPos2,...
%                     ptb_window,iPair,nPairs_current_inner,black,on_screen_instructions_line,...
%                     idx_i_real_dist)
                
                % Display the stimuli
                for iDraw = 1:nOptions
                    
                    if iDraw == 1
                        xPos_val_draw = xPos1_val;
                        yPos_val_draw = yPos1_val;
                    else
                        xPos_val_draw = xPos2_val;
                        yPos_val_draw = yPos2_val;
                    end
                    
                    switch curr_space
                        
                        case 'shade_base_space'
                            
                            % Draw auto
                            [trunk_rect_centered,base_rect_centered,...
                                shade_rect_centered,switch_rect_centered,...
                                knob_rect_centered] = draw_shade_base(...
                                SF_all,SF_body,...
                                xPos_val_draw,yPos_val_draw,op_x_loc(iDraw),op_y_loc(iDraw));
                            
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
                                xPos_val_draw,yPos_val_draw,...
                                op_x_loc(iDraw),op_y_loc(iDraw));
                            
                            Screen('DrawTextures',ptb_window,[Body_idx, Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx], [], [rBody; rBeak; rTail; rBeakTip; rTailTip]');
                            
                        case 'flower_space'
                            
                            [rStem, rStemTip, rBody, rHeart] ...
                                = drawFlower(SF_all,SF_body,...
                                xPos_val_draw,yPos_val_draw, ...
                                op_x_loc(iDraw),op_y_loc(iDraw));
                            
                            Screen('DrawTextures', ptb_window, [Stem_idx, StemTip_idx Body_idx, Heart_idx], [], [rStem; rStemTip; rBody; rHeart]');
                            
                        case 'neck_legs_space'
                            
                            [rBody, rNeck, rHead, rLegs, rFeet] ...
                                = draw_neck_legs(SF_all,SF_body,...
                                xPos_val_draw,yPos_val_draw,...
                                op_x_loc(iDraw),op_y_loc(iDraw));
                            Screen('DrawTextures', ptb_window, [Head_idx, Neck_idx Body_idx Legs_idx Feet_idx], [], [rHead; rNeck; rBody; rLegs; rFeet]');
                            
                        case 'square_circle_space'
                            
                            [rectColor,centeredRect, ovalColor, centeredOval] = ...
                                draw_square_circle(xPos_val_draw,yPos_val_draw,op_x_loc(iDraw),op_y_loc(iDraw),...
                                                   circle_x_loc(iDraw),circle_y_loc(iDraw));
                            
                            Screen('FillRect', ptb_window, rectColor, centeredRect);
                            Screen('FillOval', ptb_window, ovalColor, centeredOval);
                              
                        case 'pot_leaf_space'
                            
                            % Test
                            draw_pot_leaf(ptb_window,SF_all,SF_body,...
                                flower_idx,leaf_idx,...
                                xPos_val_draw,yPos_val_draw,op_x_loc(iDraw), op_y_loc(iDraw));
                            
                        case 'handle_circle_space'
                            
                            [propertiesMat,centeredHandle,centeredRod,...
                                rectColor,baseOval,centeredOval] = ...
                                draw_handle_circle(ptb_window,xPos_val_draw,yPos_val_draw,...
                                op_x_loc(iDraw),op_y_loc(iDraw),...
                                circle_x_loc(iDraw),circle_y_loc(iDraw));
                            
                    end
                end % iDraw
                
                % Show the mouse and keep checking its coordinates and clicks
        %             ShowCursor('Arrow');
%                 HideCursor;
 
                % Get the current position of the mouse
                [mx, my, buttons] = GetMouse(ptb_window);
                
%                 % Temp
                buttons = [1 0 0];
                pause(0.05);
                
                % Highlight the dot that mouse is closest to.
                distances = abs(scaleX - mx);
                [minDist, idx] = min(distances);

%                 % Temp record perfect responses
%                 idx = idx_i_real_dist;
                
                centeredRect = CenterRectOnPoint(baseRect_mouse,scaleX(idx),scaleY(idx));

                % Draw the rect to the screen
%                 Screen('FillRect', ptb_window, [1 0 0], centeredRect);

%                 Screen('Flip',ptb_window);

                % If a click happens, check what was clicked.
                if sum(buttons) > 0

                    response = idx;

                    mouseClickTime = GetSecs;
                    iRT = mouseClickTime - trialStartTime;

                    % Wait for mouse click release
                    while any(buttons)
                        [~, ~, buttons] = GetMouse(ptb_window);
                    end

                    % End the while loop
                    break;

                else
                    idx = NaN;
                    iRT = NaN;
                end

                % If escape pressed, exit
                [keyIsDown,~,keyCode] = KbCheck;

                if keyIsDown && keyCode(escapeKey)
                    break;
                end

            end % while

            %% Save a screenshot of the trials
            switch curr_space
                
                case 'shade_base_space'
                    
                    % GetImage call. Alter the rect argument to change the location of the screen shot
                    rect_to_save = [op_x_loc(1) - scaled_space_x(250)/2 - 5, ...
                        op_y_loc(1) - trunk_height/2 - shade_height - 5, ...
                        op_x_loc(2) + scaled_space_x(250)/2 + 5, ...
                        op_y_loc(2) + trunk_height/2 + scaled_space_y(250) + 5];
                    
                case 'pot_leaf_space'
                    
                    % GetImage call. Alter the rect argument to change the location of the screen shot
                    rect_to_save = [op_x_loc(1) - scaled_space_x(250)/2 + overlap, ...
                        op_y_loc(1) - pot_height/2 - trunk_height - flower_height - overlap, ...
                        op_x_loc(2) + scaled_space_x(250)/2 + overlap, ...
                        op_y_loc(2) + pot_height/2];

                case 'neck_legs_space'
                    
                    % GetImage call. Alter the rect argument to change the location of the screen shot
                    rect_to_save = [op_x_loc(1) - bird_width/2 + 5, ...
                                    op_y_loc(1) - bird_height/2 - head_height - scaled_space_x(250) - overlap - 3, ...
                                    op_x_loc(2) + bird_width/2 + overlap + 75, ...
                                    op_y_loc(2) + feet_height + bird_height/2 + scaled_space_x(250) - 8];
                    
                case 'beak_tail_space'
%                     
%                     % EDIT THIS FOR NEW STIMULI
%                     
                    % GetImage call. Alter the rect argument to change the location of the screen shot
                    rect_to_save = [op_x_loc(1) - bird_width/2 - beak_body_offset_x - scaled_space_x(250) + overlap - beakTip_width - 5, ...
                                    op_y_loc(1) - bird_height/2 - 5, ...
                                    op_x_loc(2) + bird_width/2 - overlap + scaled_space_y(250) + tailTip_width + 5, ...
                                    op_y_loc(2) + bird_height/2 + 5];

                case 'square_circle_space'
                    
                    % GetImage call. Alter the rect argument to change the location of the screen shot
                    rect_to_save = [op_x_loc(1) - square_diameter/2 - 5, ...
                                    op_y_loc(1) - square_diameter*coeff/2 - 5, ...
                                    op_x_loc(2) + square_diameter/2 + 5, ...
                                    circle_y_loc(2) + scaled_space_y(250)/2 + 5];             
                                
                case 'handle_circle_space'
                    
                    rect_to_save = [op_x_loc(1) - handle_width/2 - 5, ...
                        circle_y_loc(1) - scaled_space_y(250)/2 - 5,...
                        circle_x_loc(2) + scaled_space_y(250)/2 + 5,...
                        circle_y_loc(2) + scaled_space_y(250)/2 + 5];
                                        
            end % switch
            
            imageArray = Screen('GetImage', ptb_window, rect_to_save, 'backBuffer');
            
            pause(0.2);
            
            % imwrite is a Matlab function, not a PTB-3 function
            if save_data
                imwrite(imageArray, fullfile(saveLoc,saveFolder,...
                    ['pairs_' ...
                    int2str(pairs_current_inner(iPair,1)) '_' ...
                    int2str(pairs_current_inner(iPair,2)) '.png']));
            end
            
            
            %%
            if keyIsDown && keyCode(escapeKey)
                break;
            end

            % Temp
            pause(0.05);
            
            Screen('Flip',ptb_window);
            
            %% Present the ITI page

            % Set mouse position to the center
            SetMouse(xCenter, my, ptb_window);

            % Draw the scale and text
%             drawNumScale(ptb_window,scaleX,scaleY,button_size,black);
%             drawExtraStuff(pairs_current_inner(iPair,:),xPos1,yPos1,xPos2,yPos2,...
%                 ptb_window,iPair,nPairs_current_inner,black,on_screen_instructions_line,...
%                 idx_i_real_dist,idx)
            
            Screen('Flip',ptb_window);
            WaitSecs(ITI);

            KbReleaseWait;
            
        end % iPair   

end % run_dissimilarity

%% Extra functions

function drawNumScale(ptb_window,scaleX,scaleY,button_size,black)

    Screen('TextSize',ptb_window,12); 

    DrawFormattedText(ptb_window,['  lowest' '\n\ndissimilarity'],  scaleX(1) - 30, scaleY(1) - 50, black);
    DrawFormattedText(ptb_window,['  highest'  '\n\ndissimilarity'], scaleX(end) - 30, scaleY(1) - 50, black);

    % Draw all of our dots to the screen in a single line of code
    Screen('DrawDots', ptb_window, [scaleX; scaleY],...
        repmat(button_size,1,scale_max_val), black, [], 1);

    Screen('TextSize',ptb_window,15);

    % Draw the scale
    for iNum = 1:scale_max_val
        DrawFormattedText(ptb_window,int2str(iNum), scaleX(iNum) - 7, scaleY(iNum) + 35, black);
    end

end

function drawBoundaries(allCoords,x_cent,y_cent,line_width)

    % Draw the boundaries in black
    Screen('DrawLines', ptb_window, allCoords,...
        line_width, black, [x_cent y_cent]);

end

function drawExtraStuff(current_pairs,xPos1,yPos1,xPos2,yPos2,...
        ptb_window,pairDisplayed,nPairs_inner,black,on_screen_instructions_line,...
        idx_i_real_dist,idx_inner)

%     if ~exist('idx_inner','var')
%         idx_inner = 0;
%     end
%     
      Screen('TextSize',ptb_window,60);
% 
%     % Draw all text
%     DrawFormattedText(ptb_window,['Trial ' int2str(pairDisplayed) ' out of ' int2str(nPairs_inner)],20,50);
%     
      if print_corr_resp
          % Draw correct response
          DrawFormattedText(ptb_window,[int2str(idx_i_real_dist)],...
              xCenter - 120,...
              yCenter + 220);
      end
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