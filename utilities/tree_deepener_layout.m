function [XT]=tree_deepener_layout(XT,prct)

% this function makes the layout style of the hierarchy more pronounced and
% should be used after adding in one set of "verticle" edges between each
% level of the hierarchy tree represented by ClustCell

% it is assumed that the verticle distance between each unique y value in
% the tree layout XT is 1 since this is the way Matlab outputs a tree
% layout when using the "layered method"

y = XT(:,2);
yvals = unique(y);
maxy = max(yvals);

for j = 2:2:(maxy-1)
    locs = y == yvals(j);
    XT(locs,2) = yvals(j) - prct;
end