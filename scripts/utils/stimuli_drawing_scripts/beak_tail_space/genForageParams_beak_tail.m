function [SF_all,SF_body,SF_dim,rotationAngle, xPos, yPos, delta_translation, delta_rotation, dist_criterion, ...
          stand_space_x, stand_space_y, scaled_space_x, scaled_space_y, ...
          max_navig, lower_edge_idx, upper_edge_idx] ...
          = genForageParams_beak_tail
% Function to generate parameters for dials and dimentions.
% Makes sure same ones are used across experiments.

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
    
end