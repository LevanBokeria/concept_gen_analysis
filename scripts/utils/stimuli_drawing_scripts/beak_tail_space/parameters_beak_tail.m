function [beak_height, beakTip_width, beakTip_height, ...
          tail_height, tailTip_width, ...
          bird_height, bird_width, ...
          beak_body_offset_y, beak_body_offset_x, ...
          tail_body_offset, beakTip_offset, tailTip_offset, ...
          overlap] ...
          = parameters_beak_tail(SF_all,SF_body)
% This function just creates the parameters for the bird space, and passes
% it to the parent function. This makes it efficient, so that these lines
% do not have to be re-written every time.
%% Setup bird parameters

% Scaling factor. Because we'll be displaying 4 birds, each needs to be
% smaller than usual.
if (nargin < 1)
    SF_all = 1;
    SF_body = 1;
%     SF_dim  = 1;
end

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


end