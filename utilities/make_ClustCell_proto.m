function [ClustCell] = make_ClustCell_proto(Cell)

% INPUTS:
%       Cell is a row cell vriable such that each entry Cell{i} is itself a
%       row cell variable representing a clustering of the data in the
%       standard way as a row of data labels in 1:n where n is the total
%       number of data. Cell is a hierarchical system of clusterings.
%       IMPORTANTLY the order must be such that Cell{1} is the highest
%       level clustering...that is, the fewest clusters.

% the number of clusterings
num = size(Cell,2);
% initialize ClustCell
% get the number of clusters in the lowest clustering in the hierarchy
sz = size(Cell{num},2);
ClustCell = cell(sz,num+1);
ClustCell(:,1) = Cell{num}';
col = 1;
for jj = num:-1:2
    % update the column of ClustCell we will write into
    col = col + 1;
    Comms2 = Cell{jj-1};
    sz2 = size(Comms2,2);
    % initialize the storage in ClustCell
    for ii = 1:sz2
        ClustCell{ii,col} = [];
    end
    Comms1 = Cell{jj};
    sz1 = size(Comms1,2);
    % work through each cluster in Comms to locate an arbitrary point in
    % Comms2
    for ii = 1:sz1
        % get an arbitrary point
        pt = Comms1{ii};
        pt = pt(1);
        [clusters]=find_datapoint_int(Comms2,pt)';
        for xx = 1:length(clusters)
            ClustCell{clusters(xx),col} = sort([ClustCell{clusters(xx),col},ii]);
        end
    end
end

% now add in the final cell that represents the root of the hierarchy tree
sz = size(Cell{1},2);
ClustCell{1,num+1} = 1:sz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           BELOW ARE CALLED BY ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [clusters]=find_datapoint_int(CommunityCell,point)


sz=size(CommunityCell,2);
clusters=zeros(sz,1);
c=1;
for i=1:sz
    if any(ismember(point,CommunityCell{i}))
        clusters(c)=i;
        c=c+1;
    end
    
end
ind=find(clusters==0,1);
clusters(ind:end)=[];

if isempty(clusters)
    clusters=0;
end