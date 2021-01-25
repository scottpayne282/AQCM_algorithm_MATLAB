function [T,levelsizes,levels]=tree_from_cell_grid(...
    ClustCell,n)

%
% The output T is a directed network 
% which is a dendrogram of the hierarchical clustering.

% 
% the following returns a LOGICAL array of same dimensions as ClustCell
% it places a logical 1 in any location where ClustCell is nonempty
nonemptycells=~cellfun('isempty',ClustCell);
%
levelsizes=sum(nonemptycells,1); % returns # clusters in each level as row
%
levels=size(nonemptycells,2); % returns the number of levels
%
nodesum=sum(levelsizes)+n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=sparse(nodesum,nodesum);

%
currentnode=nodesum;
for j=levels:-1:1;
    if j>1;
        % in this step the outer loop on j has just "stepped down" by one
        % column (level) in the ClustCell data.
        % Since the inner loop on i will cout backwards as well, this means
        % that the "currentnode" we are working on is represented by the
        % bottom block in level j of ClustCell.
        % The following code calculates the row array "prevlev" which has
        % the vertex labels (for the graph T) of the clusters represented
        % in ClustCell by level j-1.
        % In the "i loop" below we will fill in the edges of T between
        % level j and level j-1.
        backslide=sum(levelsizes([j-1 j]))-1;
        prevlev=((currentnode-backslide):(currentnode-levelsizes(j)));
    else
        prevlev=(1:n);
    end    
    for i=levelsizes(j):-1:1;
        % We must:
        %       Access ClustCell{i,j}
        %       Replace the labels in that cell with the ones from 
        %       the "prevlev" list we made above: Just use
        %       prevlev like a function with input as the list we got
        %       out of the ClustCell{i,j}
        %       
        %       The edges we want to add to our list are:
        %       currentnode->each new label we just got from the above
        %       steps.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        labels=ClustCell{i,j};
        sz=size(labels,2);
        %
        % now we may store the edges to the graph T
        T(currentnode,prevlev(labels))=ones(1,sz);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currentnode=currentnode-1;
    end
end
% 

