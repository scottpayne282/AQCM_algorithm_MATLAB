function [T,Tlabels,CommunityCell]=community_detectMINavMemSparse(...
    ClustCell,InfoCell,n)

%-------------------------------------------------------------------------
[T,E,root,sizelist,~,~]=dendrogramerMAXtreeMinav(...
    ClustCell,InfoCell,n);
% E is the 3 column edge list of the maximum spanning tree of T
% sizelist will be used to "trim" the branches of the tree represented by
% E which correspond to single vertex clusters.
sizelist=sizelist(n+1:end);
sizelist(end)=[];
%-------------------------------------------------------------------------
% in here we "trim the leaves" of E
E=E(n+1:end,:);
%-------------------------------------------------------------------------
% NOTE that the above lines adjusting the dimensions of "sizelist" and "E" 
% allow the two arrays to correspond, meaning that each row of E is an edge
% and the column of sizelist of the same index represents the size of the
% cluster which is the head vertex of that edge.
%-------------------------------------------------------------------------
% in here we trim the single vertex cluster branches.
E(sizelist==1,:)=[];
%-------------------------------------------------------------------------
[Tlabels]=minavgcutInside(E,root);
sz=size(Tlabels,2);
Tlabels=sort(Tlabels);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% the following code stores (to a cell) the numeric labels (of the vertices
% of the original graph) of the clusters 
% 
[CommunityCell]=leaves_from_tree_nodes_inside(T,n,Tlabels);
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               BELOW ARE CALLED BY THE ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [T,E,root,sizelist,levelsizes,levels]=dendrogramerMAXtreeMinav(...
    ClustCell,InfoCell,n)

%
% The output T is a directed network 
% which is a dendrogram of the hierarchical clustering.
% E is the maximum spanning tree stored as an edge list with labels
% levelsizes, levels are tools used WITHIN dendrogramer to built the tree T
% out of the data ClustCell. The reason we output those is that we use them
% in community_detect to build the coordinates X in order to plot a diagram
% of T
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
root=nodesum;
T=sparse(nodesum,nodesum);
sizelist=ones(1,nodesum);
neglocs = [];
zerolocs = [];
nanlocs = [];
inflocs = [];
minT = Inf;
maxT = 0;
sumT = 0;
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
        
        if j>1;
            % create a weight column corresponding to arcs
            % goto the earlier "level j-1" and get the desired data.
            % It is VERY important to note that the following command must be run by
            % itself BEFORE attempting to treat the object as an array....in other
            % words we store to "data" first so that we will be able to use array
            % syntax to extract info
            data=cell2mat(InfoCell(labels,j-1));
            fatherden=InfoCell{i,j}(1); % This works...I checked it
            %
            % We can count by 2's to get the data separated...recall that in any cell
            % in InfoCell the array stored there is a 2x1 where the top entry is a
            % density of the cluster represented by that cell and the bottom entry is
            % the size of that cluster.
            sonden=(data(1:2:(2*sz)))';
            sonverts=(data(2:2:(2*sz)))';
            % We use element by element operations to output all of the
            % calculations simultaneously as a row
            %
            arcwghts=sonverts./((sonden.^2)-(fatherden^2));
            % here we track indices of some possible special cases
            if length(arcwghts)==1
                arcwghts = Inf;
            end
            parent = currentnode;
            children = prevlev(labels);
            % case of negative
            locs = arcwghts<0;
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                neglocs = [neglocs; [repmat( parent,[numt 1] ) temp ] ];
            end
            % case of zero
            locs = arcwghts==0;
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                zerolocs = [zerolocs; [repmat( parent,[numt 1] ) temp ] ];
            end
            % case of NaN
            locs = isnan(arcwghts);
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                nanlocs = [nanlocs; [repmat( parent,[numt 1] ) temp ] ];
            end
            % case of Inf
            locs = arcwghts==Inf;
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                inflocs = [inflocs; [repmat( parent,[numt 1] ) temp ] ];
            end
            % now track min's and max's and sum
            minim = min(arcwghts(arcwghts>0));
            if ~isempty(minim)
                minT = min(minT,minim);
            end
            maxim = max(arcwghts(arcwghts~=Inf));
            if ~isempty(maxim)
                maxT = max(maxT,maxim);
            end
            sumim = sum(arcwghts(arcwghts>0 & arcwghts~=Inf));
            if ~isempty(sumim)
                sumT = sumT + sumim;
            end
            sizelist(prevlev(labels))=sonverts;
        else
            arcwghts=Inf(1,sz); % the edges between level 1 clusters and single 
            % vertices have flag value "-1" so we can catch them below and
            % switch them to high values
            parent = currentnode;
            children = prevlev(labels);
            % case of Inf
            locs = arcwghts==Inf;
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                inflocs = [inflocs; [repmat( parent,[numt 1] ) temp ] ];
            end
        end
        % now we may store the edges to the graph T
        T(currentnode,prevlev(labels))=arcwghts;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currentnode=currentnode-1;
    end
end
% 

sizelist(end)=n;
% below we will do the switching of flagged values
%-------------------------------------------------------------------------
if ~isempty(inflocs)
    for iii = 1:size(inflocs,1)
        T(inflocs(iii,1),inflocs(iii,2)) = 2*sumT;
    end
end
%
if ~isempty(neglocs)
    for iii = 1:size(neglocs,1)
        T(neglocs(iii,1),neglocs(iii,2)) = 2*sumT;
    end
end
%
if ~isempty(zerolocs)
    for iii = 1:size(zerolocs,1)
        T(zerolocs(iii,1),zerolocs(iii,2)) = (1/n)*minT;
    end
end
%
if ~isempty(nanlocs)
    for iii = 1:size(nanlocs,1)
        T(nanlocs(iii,1),nanlocs(iii,2)) = (1/n)*minT;
    end
end
%-------------------------------------------------------------------------
% Now that we have complete information about the graph T we may store the
% maximum spanning tree to the edgelist E
E=zeros(nodesum-1,3);
for j=1:(nodesum-1)
    maxj=full(max(T(:,j)));
    E(j,:)=[maxj find(T(:,j)==maxj,1) j];
    % observe that the 3rd column is the vertex label of the head of the
    % edge
end
%-------------------------------------------------------------------------
% 
if ~isempty(inflocs)
    for iii = 1:size(inflocs,1)
        T(inflocs(iii,1),inflocs(iii,2)) = maxT;
    end
end
%
if ~isempty(neglocs)
    for iii = 1:size(neglocs,1)
        T(neglocs(iii,1),neglocs(iii,2)) = maxT;
    end
end
% the reason for this is to fix the scale of the edge 
% weights of T(L). We had set the very low for the purpose of our community
% selection method, however when we display the dendrogram T we want to be
% able to visually see the edge weights on a color scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [labels]=minavgcutInside(E,root)
%
% INPUT: root is the vertex label of the root vertex of the tree
% E is an edge list (n-1x3 array) where n is the number of vertices
% in the graph (tree) represented by E.
% column 1 is the edge weight
% column 2 is the tail vertex of the edge
% column 3 is the head vertex of the edge
% The edge label is assumed to be the HEAD vertex of the edge
% * NOTE that the smallest label NEED NOT be 1 since we may be dealing with
% a subtree of a larger tree or some graph where the vertex labels do not
% start with 1
%
% OUTPUT: "labels" is a list of the vertex labels of the head vertices of
% the edges which are a root separating cut of maximum average weight
%
% First we must add some information to E, so we ADD a column in frontof E
% The new first column will be the lambda values used in the algorithm.
sz=size(E,1);
E=[zeros(sz,1) E];
%-------------------------------------------------------------------------
% now we will calculate the lambda values
for e=1:sz
    if any(E(:,3)==E(e,4))
        E(e,1)=(sum(E(E(:,3)==E(e,4),2))-(E(e,2)))/...
            (length(find(E(:,3)==E(e,4)))-1);
    else
        E(e,1)=Inf;
    end
end
E(isnan(E(:,1)),1)=-Inf;
%-------------------------------------------------------------------------
% now we sort the edgelist (in highest to lowest order) E by the values 
% in column 1
E=sortrows(E);
%-------------------------------------------------------------------------
% calculate the average weight of the cut "located at" the root
alph0=mean(E(E(:,3)==root,2));
%-------------------------------------------------------------------------
% now we perform the contraction based algorithm
while E(1,1)<alph0
    if E(1,3)==root
        % here we will contract the edge and update alph0
        % all children edges of the contracted edge have their tail reassigned
        % to the tail vertex of that contracted edge
        E(E(:,3)==E(1,4),3)=E(1,3);
        % delete the contracted edge
        E(1,:)=[];
        alph0=mean(E(E(:,3)==root,2));
    else
        % here we will contract the edge and update the lambda value of the
        % edge that was the "in edge" of our edge. Then we will need to
        % "re-file" the updated edge so that E remains sorted by lambda
        % order
        % re-assign children (same code as above)
        E(E(:,3)==E(1,4),3)=E(1,3);
        % store the tail vertex of the edge we are trying to contract
        tail=E(1,3);
        % delete the contracted edge
        E(1,:)=[];
        % recalculate the lambda value for the "in edge" to the one we just
        % contracted
        e=find(E(:,4)==tail,1);
        E(e,1)=(sum(E(E(:,3)==E(e,4),2))-(E(e,2)))/...
            (length(find(E(:,3)==E(e,4)))-1);
        % now re-locate edge "e" so that the list E remains lambda sorted
        if E(1,1)>E(e,1)
            E=[E(e,:);E];
            E(e+1,:)=[];
        else
            erow=E(e,:);
            E(e,:)=[];
            loc=find(E(:,1)>erow(1),1);
            if isempty(loc)
                E=[E;erow];
            else
                E=[E((1:loc-1),:);erow;E(loc:end,:)];
            end
            %
        end
    end
end
% now we record the head vertices of the edges whose tail is at the root
labels=(E(E(:,3)==root,4))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Comms]=leaves_from_tree_nodes_inside(T,n,roots)

% assume indexes 1:n is the bottom row (bottom level) of the rooted
% hierarchy tree T (all arrows point down)

% the input variable "roots" is a row array of numeric labels of nodes of
% the tree T for which we wish to find the leaves

% the stack variable stk stores items by row where the first entry is the
% node number of a node and the second is it's level inthe hierarchy tree

num=size(T,1);
sz = size(roots,2);
Comms = cell(1,sz);

for ii=1:sz
    % initialize storage for the leaves
    leaves = false(1,n);
    % initialize a stack
    stk=zeros(num,1);
    hgt=1; % current stack height
    % initialize the stk with the root as the only item
    stk(1)=roots(ii);
    % initialize a variable to keep track of previous children nodes
    prev_chld = [];
    while hgt>0
        % get a node from the top of the stack
        node=stk(hgt);
        stk(hgt)=0;
        hgt=hgt-1;
        % get the children of that node
        chld=find(T(node(1),:));
        nchld=size(chld,2);
        if nchld>0
            % record the current children
            prev_chld = chld;
            % push the children to the stack
            stk((hgt+1):(hgt+nchld))=chld';
            % update the stack height
            hgt=hgt+nchld;
        else
            leaves(prev_chld) = true;
        end
    end
    Comms{ii} = find(leaves);
end
% 







