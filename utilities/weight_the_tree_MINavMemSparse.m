function [T,E,root,sizelist]=weight_the_tree_MINavMemSparse(T,S)

n = size(S,1);

% initialize levjmin1 as the leaves
levjmin1 = find(~any(T,2))';
levj = find(any(T(:,levjmin1),2))';

nodesum = size(T,1);
root = nodesum;
% this variable is used to keep track of the sizes of the clusters that the
% tree nodes represent
sizelist=ones(1,nodesum);

%--------------------------------------------------------------------------
neglocs = [];
zerolocs = [];
nanlocs = [];
inflocs = [];
minT = Inf;
maxT = 0;
sumT = 0;
%--------------------------------------------------------------------------

szlevj = length(levj);
% initialize some necessary storage
ClustCell = cell(szlevj,ceil(0.5*nodesum));
InfoCell = cell(szlevj,ceil(0.5*nodesum));

% handle the leaf edges first...they will each have weight Inf
level = 1;
for i = 1:szlevj
    currnode = levj(i);
    labels = find(T(currnode,:));
    % record the labels of the data in the cluster
    ClustCell{i,level} = labels;
    % get the size of the cluster
    sizec = length(labels);
    % calculate the density of the cluster using the code identical to that
    % used in the QCMP software
    if sizec>1
        temp=nchoosek(labels',2);
        lininds=sub2ind([n,n],temp(:,1),temp(:,2))';
        InfoCell{i,level} = [2*sum(S(lininds))/(sizec*(sizec-1));sizec];
    else
        InfoCell{i,level} = [Inf;sizec];
    end
    T(currnode,labels) = Inf;
            parent = currnode;
            children = labels;
            % case of Inf
            locs = true(1,sizec);
            if any(locs)
                temp = children(locs)';
                numt = length(temp);
                inflocs = [inflocs; [repmat( parent,[numt 1] ) temp ] ];
            end
end

stillworking = true;

while stillworking
    level = level + 1;
    levjmin1 = levj;
    levj = find(any(T(:,levjmin1),2))';
    if length(levj)==1
        stillworking = false;
    end
    szlevj = length(levj);
    for i = 1:szlevj
        currnode = levj(i);
        labels = find(T(currnode,:));
        % record the labels of the data in the cluster
        currcluster = [];
        for ii = 1:length(labels)
            currcluster = union(currcluster,ClustCell{ levjmin1==labels(ii)  ,level-1});
        end
        if length(currcluster)>size(currcluster,2)
            currcluster = currcluster';
        end
        ClustCell{i,level} = currcluster;
        % get the size of the cluster
        sizec = length(currcluster);
        % calculate the density of the cluster using the code identical to that
        % used in the QCMP software
        if sizec>1
            temp=nchoosek(currcluster',2);
            lininds=sub2ind([n,n],temp(:,1),temp(:,2))';
            InfoCell{i,level} = [2*sum(S(lininds))/(sizec*(sizec-1));sizec];
        else
            InfoCell{i,level} = [Inf;sizec];
        end
        % calculate and write edgeweights to T using the same technique in
        % standard T building methods of QCMP software package
        locs = ismember(levjmin1,labels);
        data=cell2mat(InfoCell(locs,level-1));
        fatherden=InfoCell{i,level}(1);
        sz = length(labels);
        sonden=(data(1:2:(2*sz)))';
        sonverts=(data(2:2:(2*sz)))';
        arcwghts=sonverts./((sonden.^2)-(fatherden^2));
        if length(arcwghts)==1
            arcwghts = Inf;
        end
        %-------------------------------------------------
        % here we track indices of some possible special cases
            parent = currnode;
            children = labels;
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
        %-------------------------------------------------
        sizelist( labels )=sonverts;
        T(currnode, labels )=arcwghts;
    end
end

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
    maxj=max(T(:,j));
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








