function [Beak, BeakTip, Tail, TailTip, Body, Bar_Beak, Bar_Tail, ...
          Beak_idx, Tail_idx, BeakTip_idx, TailTip_idx, Body_idx, ...
          Bar_Beak_idx, Bar_Tail_idx] ...
          = load_beak_tail_parts(home,window)
% This function just loads the bird parts and makes them into textures. Then passes
% it to the parent function. This makes it efficient, so that these lines
% do not have to be re-written every time.
%% Load the bird parts

% load the bits of the bird
Beak    = imread(fullfile(home,'inputs','beak_tail_space','birdparts','beak.png'),'BackgroundColor',[1,1,1]);    % Beak(:,:,4) = Beak(:,:,1);
BeakTip = imread(fullfile(home,'inputs','beak_tail_space','birdparts','beak_tip.png'),'BackgroundColor',[1,1,1]); % BeakTip(:,:,4) = BeakTip(:,:,1);
Tail    = imread(fullfile(home,'inputs','beak_tail_space','birdparts','tail.png'),'BackgroundColor',[1,1,1]);    % Tail(:,:,4) = Tail(:,:,1);
TailTip = imread(fullfile(home,'inputs','beak_tail_space','birdparts','tail_tip.png'),'BackgroundColor',[1,1,1]); % TailTip(:,:,4) = TailTip(:,:,1);
Body    = imread(fullfile(home,'inputs','beak_tail_space','birdparts','body.png'),'BackgroundColor',[1,1,1]);    % Body(:,:,4) = Body(:,:,1);

% bars for choosing the neck:legs ratio
Bar_Beak = imread(fullfile('inputs','beak_tail_space','birdparts','Bar_controller.png')); Bar_Beak(:,:,4) = 255;                                     
Bar_Tail = imread(fullfile('inputs','beak_tail_space','birdparts','Bar_controller.png')); Bar_Tail(:,:,4) = 255;                                       
 
%% Make them into textures
Beak_idx    = Screen('MakeTexture', window, Beak);
Tail_idx    = Screen('MakeTexture', window, Tail);
BeakTip_idx = Screen('MakeTexture', window, BeakTip);
TailTip_idx = Screen('MakeTexture', window, TailTip);
Body_idx    = Screen('MakeTexture', window, Body);

Bar_Beak_idx = Screen('MakeTexture', window, Bar_Beak);
Bar_Tail_idx = Screen('MakeTexture', window, Bar_Tail); 


end