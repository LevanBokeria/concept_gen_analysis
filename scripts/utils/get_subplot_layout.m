% Functions that creates a logical arrangement of nr_plots in one figure
% Last edited: 14-11-2014 by Erik Meijs

function [nr_rows, nr_columns] = get_subplot_layout(nr_plots)
% Input:
%   nr_plots        =   number of plots you want to fit into one figure
%
% Outputs:
%   nr_rows         =   number of rows of subplots in the figure
%   nr_columns      =   number of columns of subplots in the figure


%% CALCULATE VALUES
root_nr_plots   = sqrt(nr_plots);               % Get square root (approximation of number of plots)

nr_rows         = floor(root_nr_plots);         % Get number of rows (lower amount of square root)
nr_columns      = ceil(nr_plots/nr_rows);       % Get number of columns based on total amount and number of rows