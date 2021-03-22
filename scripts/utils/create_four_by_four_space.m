function [G,coordinate_matrix, unit_dist] = create_four_by_four_space

% Make a figure to plot triangles
D = zeros(16,16);

connected_nodes = [1,2;2,3;3,4;5,6;6,7;7,8;9,10;10,11;11,12;13,14;14,15;...
    15,16;1,5;5,9;9,13;2,6;6,10;10,14;3,7;7,11;11,15;4,8;8,12;12,16];

for iCon = 1:size(connected_nodes,1)
    D(connected_nodes(iCon,1),connected_nodes(iCon,2)) = 1;
    D(connected_nodes(iCon,2),connected_nodes(iCon,1)) = 1;
end

G = graph(D,sprintfc('%d',1:16));

coordinate_matrix = [100 100 100 100 ...
                     150 150 150 150 ...
                     200 200 200 200 ...
                     250 250 250 250;...
                     repmat([100 150 200 250],1,4)]';

unit_dist = 50;                 
                 
% h = plot(G);
% h.XData = coordinate_matrix(:,1);
% h.YData = coordinate_matrix(:,2);

end