function find_good_companions_for_triangles(source_triangles)

%% Description

% This function will find good matching triangles for a given triangle, in
% our 4x4 stimulus space. It accepts multiple conditions as flags, such as
% whether to exclude certain nodes, whether the companion triangle can
% share dimensions with the original once, etc etc.

% The function takes a matrix that specifies nodes of all the triangles for
% which companion triangles must be found, and loops through them producing
% corresponding list of companion triangles.

%% General setup
dbstop if error;
% close all;

% Check you're in the right directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/unified_space/');
end

addpath(genpath(home));

if (nargin < 1) 
    source_triangles = [14,9,4];
end


make_triangle_plots = 1;
remove_node_labels  = 0;
save_figures        = 0;

[G,coordinate_matrix,unit_dist] = create_four_by_four_space;

all_nodes = 1:16;

all_node_combs = combnk(all_nodes,3);

%% %%%%%%% Conditions for the companion triangle %%%%%%%%

% Which nodes to exclude
exclude_nodes_base = [1,2,5,6,16];
% exclude_nodes_base = [1,16];

% Which nodes to necessarily have?
include_nodes = [];

% Minimum city block distance between elements
min_cb_dist = 2;

% Can two toys share a dimension?
can_share_dimension = 0;

% Must two toys share a dimension?
must_share_dimension = 0;

% Must it have at least one target in the center?
must_contain_center = 1;

% Must it have at least one target on the edge?
must_contain_edge = 0;

% Must it NOT have any targets in the center?
must_avoid_center = 0;

%% Condition relating to the source triangle

% Can they share a node?
must_avoid_same_nodes_with_source = 1;

% Can they share a dimension?
must_avoid_same_dimensions_with_source = 0;

% If you take the nodes that are opposite of the source triangle nodes 
% across the 45 degree line, can those nodes be shared?
must_avoid_mirrored_nodes_with_source = 1;

%% Start looping through the given list of triangles

for iTr = 1:size(source_triangles,1)
    
    exclude_nodes = exclude_nodes_base;
    
    source_tr_nodes = source_triangles(iTr,:);
    
    % depending on the source triangle conditions, update the companion
    % triangle conditions
    if must_avoid_same_nodes_with_source
        exclude_nodes = [exclude_nodes, source_tr_nodes];
    end
    
    if must_avoid_mirrored_nodes_with_source
       
        % Find the mirror nodes
        source_tr_coords = coordinate_matrix(source_tr_nodes,:);
        
        % 1. Flip the coordinates of the current nodes:
        source_tr_coords_mirrored = fliplr(source_tr_coords);
        
        % 2. Find the node indices of the flipped coordinates
        [~, source_tr_nodes_mirrored] = ismember(...
            source_tr_coords_mirrored,coordinate_matrix,'rows');
        source_tr_nodes_mirrored = source_tr_nodes_mirrored';
        
        % 3. Add these to the excluded nodes
        exclude_nodes = [exclude_nodes, source_tr_nodes_mirrored];
    end
    
    %% %%%%%%%%%% Filter out the combinations %%%%%%%%%%%%%%%
    
    % Vector to signify which combs to filter
    filter_comb_logical = zeros(size(all_node_combs,1),1);
    
    for iComb = 1:size(all_node_combs,1)
        
        curr_comb = all_node_combs(iComb,:);
        
        % For every pair in the triplet, check the conditions
        all_pairs = combnk(curr_comb,2);
        
        % Pregen a logical variables for conditions
        some_pairs_share_dim = 0;
        comb_contains_center = 0;
        comb_contains_edge   = 0;
        comb_contains_included_node = 0;
        comb_shared_dim_with_source = 0;
        
        for iNode = 1:size(curr_comb,2)
            
            curr_coords = coordinate_matrix(curr_comb(iNode),:);
            
            % Is it an excluded node
            if ~isempty(find(ismember(curr_comb,exclude_nodes)))
                filter_comb_logical(iComb) = 1;
            end
            
            % Is it an included node?
            if ~isempty(find(ismember(curr_comb,include_nodes)))
                comb_contains_included_node = 1;
            end
            
            % Is it a central node?
            if isempty(find(curr_coords == 100)) && ...
                    isempty(find(curr_coords == 250))
                
                comb_contains_center = 1;
            end
            
            % Is it an edge node?
            if ~isempty(find(curr_coords == 100)) || ...
                    ~isempty(find(curr_coords == 250))
                
                comb_contains_edge = 1;
            end
        end % for iNode
        
        if ~isempty(include_nodes)
            if ~comb_contains_included_node
                filter_comb_logical(iComb) = 1;
            end
        end
        
        if must_contain_center
            if ~comb_contains_center
                filter_comb_logical(iComb) = 1;
            end
        end
        if must_contain_edge
            if ~comb_contains_edge
                filter_comb_logical(iComb) = 1;
            end
        end
        if must_avoid_center
            if comb_contains_center
                filter_comb_logical(iComb) = 1;
            end
        end
        
        %% Check the pairs for distance and sharing dimensions 
        for iPair = 1:size(all_pairs,1)
            
            curr_pair = all_pairs(iPair,:);
            
            curr_coords = coordinate_matrix(curr_pair,:);
            
            % If distance between the pairs if less than minimum, exlcude
            curr_dist = pdist(curr_coords,'cityblock');
            
            if curr_dist < min_cb_dist*unit_dist
                filter_comb_logical(iComb) = 1;
                %            break;
            end
            
            % Do the pairs share dimensions?
            component_distances = curr_coords(1,:) - curr_coords(2,:);
            
            if find(component_distances == 0)
                some_pairs_share_dim = 1;
            end
        end % iPair
        
        % If must share dimensions, and none of the pairs do:
        if must_share_dimension
            if ~some_pairs_share_dim
                filter_comb_logical(iComb) = 1;
            end
        end
        
        % If the must NOT share dimensions:
        if ~can_share_dimension
            if some_pairs_share_dim
                filter_comb_logical(iComb) = 1;
            end
        end
        
        %% Shares dimensions with the source triangle?
        
        % Find all combinations of nodes between source and companion
        [X,Y] = meshgrid(source_tr_nodes,curr_comb);
        all_source_comp_pairs = [X(:) Y(:)];
        
        for iPair = 1:size(all_source_comp_pairs,1)
            
            curr_pair = all_source_comp_pairs(iPair,:);
            
            curr_coords = coordinate_matrix(curr_pair,:);
                  
            % Do the pairs share dimensions?
            component_distances = curr_coords(1,:) - curr_coords(2,:);
            
            if find(component_distances == 0)
                comb_shared_dim_with_source = 1;
            end
        end % iPair
        
        if must_avoid_same_dimensions_with_source
            if comb_shared_dim_with_source
                filter_comb_logical(iComb) = 1;    
            end
        end
    end % iComb
    
    % Filter the combinations
    filtered_combs = all_node_combs;
    
    filtered_combs(boolean(filter_comb_logical'),:) = [];
    
    figure
%     set(gcf,'Position',[0.0010    0.0410    1.2800    0.6073]*1000);
    
%     set(gcf,'Position',[0.0957    0.0603    1.1720    0.6897]*1000);    
    [nrow,ncol] = get_subplot_layout(size(filtered_combs,1));
    
    if ncol == 2
        set(gcf,'Position',[ 444.0000   60.3000  823.7000  759.7000]);
    elseif ncol == 3
        set(gcf,'Position',[0.1500    0.0603    1.1177    0.6897]*1000);
    end
%     nrow = size(filtered_combs,1);
%     ncol = 1;
    %% Plot the triangles
    if make_triangle_plots
        
        source_tr_node_string = sprintf('%.0f_',source_tr_nodes);
        source_tr_node_string = source_tr_node_string(1:end-1);
        
        for iComb = 1:size(filtered_combs,1)
            
            % First, plot the source triangle
            subplot(nrow,ncol,iComb)

            source_coords = coordinate_matrix(source_tr_nodes,:);
            
            plot_triangle(source_coords,source_tr_nodes,G,...
                coordinate_matrix,remove_node_labels,'b','r')

            
            % Now the companion
            curr_comb = filtered_combs(iComb,:);
            
            % Draw this specific triangle
            iCoord = coordinate_matrix(curr_comb,:);
            
%             subplot(nrow,ncol,iComb*2)
            hold on
            plot_triangle(iCoord,curr_comb,G,...
                coordinate_matrix,remove_node_labels,'g','g')
            title(['#' int2str(iComb)]);
        end % iComb
        
        sgtitle(['Source tr nodes: ' ...
            source_tr_node_string],'Interpreter','none');
    end
    
    name_to_save = ['source_tr_' source_tr_node_string '_companions'];
    
%     name_to_save = ['min_dist ' int2str(min_cb_dist)...
%         '; allowed_share_dimension ' int2str(can_share_dimension)...
%         '; must share dim ' int2str(must_share_dimension)...
%         '; must center ' int2str(must_contain_center)...
%         '; must edge ' int2str(must_contain_edge)...
%         '; excl_nodes ' int2str(~isempty(exclude_nodes))];
   
    % Save
    if save_figures
        saveas(gcf,fullfile(home,'results','analysis','sandbox',...
            [name_to_save '.png']))
    end
    
end % for each source triangle

%% Subfunctions



end