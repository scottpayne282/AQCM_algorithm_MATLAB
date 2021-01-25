function [NewCell]=tree_deepener(ClustCell,depth)

% this function adds "edges" between levels in the tree described by the
% cell variable ClustCell. The purpose is to make the tree easier to
% visually analyze in situations when the tree is large and complex.

levels = size(ClustCell,2);

sz1 = size(ClustCell,1);

newsize = (depth+1)*(levels-1)+1;
NewCell = cell(sz1,newsize);

szgrid = cellfun('length',ClustCell);

for j = 1:levels-1
    % store the clustering info for level j
    newj = depth*(j-1) + j;
    NewCell(:,newj) = ClustCell(:,j);
    % calculate the number of clusters represented by the column of
    % ClustCell
    numclust = find(szgrid(:,j)>0,1,'last');
    for jj = 1:depth
        for iii = 1:numclust
            NewCell{iii,newj+jj} = iii;
        end
    end
end

j = levels;
% store the clustering info for level j
newj = depth*(j-1) + j;
NewCell(:,newj) = ClustCell(:,j);













