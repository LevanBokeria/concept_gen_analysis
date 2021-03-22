function plot_triangle(iCoord,nodes,G,coordinate_matrix,...
    remove_node_labels,tr_color,node_color)

h = plot(G);
h.XData = coordinate_matrix(:,1);
h.YData = coordinate_matrix(:,2);
h.NodeFontSize = 10;
h.MarkerSize = 1;
patch('Faces',[1,2,3],'Vertices',iCoord,...
    'FaceColor', tr_color, ...
    'FaceAlpha', 0.2);
highlight(h,nodes,'NodeColor',node_color,'MarkerSize',3)

if remove_node_labels
    labelnode(h,1:16,'');
end

end