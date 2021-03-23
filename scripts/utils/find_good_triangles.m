%% Visualize all possible triangle organizations of toys in our 4x4 space
clear; clc;
close all

dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./concept_gen_analysis/');
end

addpath(genpath(home));

make_triangle_plots      = 1;
remove_node_labels       = 0;
save_figures             = 0;
find_companion_triangles = 1;

[G,coordinate_matrix,unit_dist] = create_four_by_four_space;

all_nodes = 1:16;

all_node_combs = combnk(all_nodes,3);

%% %%%%%%% Conditions %%%%%%%%

% Which nodes to exclude
% exclude_nodes = [1,4,13,16];
exclude_nodes = [1,2,6,5,16,4,13];

% Which nodes to necessarily have?
include_nodes = [];

% Minimum city block distance between elements
min_cb_dist = 2;

% Can two toys share a dimension?
can_share_dimension = 1;

% Must two toys share a dimension?
must_share_dimension = 1;

% Must it have at least one target in the center?
must_contain_center = 1;

% Must it have at least one target on the edge?
must_contain_edge = 0;

% Must it NOT have any targets in the center?
must_avoid_center = 0;

% Among the combinations that satisfy all the above, should we exclude
% mirrored ones? So if one tr is a mirror of another across the 45 degree
% line, shall we exclude that?
exclude_mirror_triangles = 1;

%% If we will also search for comparison triangles:
% Arguments to pass that function, specifying relationship between the
% source and companion triangles.

% Can they share a node?
must_avoid_same_nodes_with_source = 1;

% Can they share a dimension?
must_avoid_same_dimensions_with_source = 0;

% If you take the nodes that are opposite of the source triangle nodes
% across the 45 degree line, can those nodes be shared?
must_avoid_mirrored_nodes_with_source = 1;
%% %%%%%%%%%% Filter out the combinations %%%%%%%%%%%%%%%

% Vector to signify which combs to filter
filter_comb_logical = zeros(size(all_node_combs,1),1);

for iComb = 1:size(all_node_combs,1)
    
    curr_comb = all_node_combs(iComb,:);
    
    % For every pair in the triplet, check the conditions
    all_pairs = combnk(curr_comb,2);
    
    % Pregen a logical recording if any pairs share dimension
    some_pairs_share_dim = 0;
    comb_contains_center = 0;
    comb_contains_edge   = 0;    
    comb_contains_included_node = 0;
    
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
    
    for iPair = 1:size(all_pairs,1)
        
        curr_pair = all_pairs(iPair,:);
        
        curr_coords = coordinate_matrix(curr_pair,:);

        %% If distance between the pairs if less than minimum, exlcude
        curr_dist = pdist(curr_coords,'cityblock');
        
        if curr_dist < min_cb_dist*unit_dist
           filter_comb_logical(iComb) = 1;
%            break;
        end

        %% Do the pairs share dimensions?
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
    
end % iComb

% Filter the combinations 
filtered_combs = all_node_combs;

filtered_combs(boolean(filter_comb_logical'),:) = [];

%% Exclude mirrored versions?
if exclude_mirror_triangles
    
    filtered_combs_no_mirror = filtered_combs;
    
    cleaning_done = 0;
    ctr = 1;
    
    while ~cleaning_done 
            
        iComb = filtered_combs_no_mirror(ctr,:)
        
        iComb_coords = coordinate_matrix(iComb,:);
       
        % 1. Flip the coordinates of the current nodes:
        iComb_coords_mirrored = fliplr(iComb_coords);
        
        % 2. Find the node indices of the flipped coordinates
        [~, iComb_nodes_mirrored] = ismember(...
            iComb_coords_mirrored,coordinate_matrix,'rows');
        iComb_nodes_mirrored = iComb_nodes_mirrored'
        
        % 3. Are these part of the filtered_combs? 
        
        % 3.1 Get all the combinations of the nodes
        mirror_node_perms = perms(iComb_nodes_mirrored);
        
        [~,mirror_idx] = ismember(mirror_node_perms,...
                                  filtered_combs_no_mirror,'rows')
        
        assert(nnz(mirror_idx) < 2);
        
        idx_to_remove = mirror_idx(mirror_idx ~= 0);
        
        % 4. If it is, delete the matching row
        if ~isempty(idx_to_remove)
            
            % Make sure that this mirror is not identical with the
            % original. Sometimes, the original had itself as the mirror!
            if ~isequal(filtered_combs_no_mirror(idx_to_remove,:), ...
                    iComb)

                filtered_combs_no_mirror(idx_to_remove,:) = [];
            end
        end
        
        
        if ctr == size(filtered_combs_no_mirror,1)
            cleaning_done = 1;
        end
        
        ctr = ctr + 1;
    
    end
    
    filtered_combs = filtered_combs_no_mirror;
end % if exclude mirror triangles



figure
set(gcf,'Position',[ 1          41        1280         907]);

[nrow,ncol] = get_subplot_layout(size(filtered_combs,1));

%% Plot each triangle
if make_triangle_plots
    for iComb = 1:size(filtered_combs,1)
        
        curr_comb = filtered_combs(iComb,:);
        
        % Draw this specific triangle
        iCoord = coordinate_matrix(curr_comb,:);
        
        subplot(nrow,ncol,iComb)
        
        plot_triangle(iCoord,curr_comb,G,coordinate_matrix,...
            remove_node_labels,'b','r')
        
        title(['# ' int2str(iComb)]);
    end % iComb
end

% name_to_save = ['min_dist ' int2str(min_cb_dist)...
%          '; allowed_share_dimension ' int2str(can_share_dimension)...
%          '; must share dim ' int2str(must_share_dimension)...         
%          '; must center ' int2str(must_contain_center)...
%          '; must edge ' int2str(must_contain_edge)...
%          '; excl_nodes ' int2str(~isempty(exclude_nodes))];

% Save
if save_figures
    name_to_save = input('Whats the name?\n');
    
    % Title
    sgtitle([name_to_save '; Total triangles = ' int2str(size(filtered_combs,1))],...
        'Interpreter','none')
    
    
    saveas(gcf,fullfile(home,'results','analysis','sandbox',...
        [name_to_save '.png']))
end

%% Find companion triangles for each? 

if find_companion_triangles
    
    find_good_companions_for_triangles(source_triangles,...
        exclude_nodes_base,include_nodes,min_cb_dist,can_share_dimension,...
        must_share_dimension,must_contain_center,must_contain_edge,...
        must_avoid_center,...
        must_avoid_same_nodes_with_source,...
        must_avoid_same_dimensions_with_source,...
        must_avoid_mirrored_nodes_with_source)
end
