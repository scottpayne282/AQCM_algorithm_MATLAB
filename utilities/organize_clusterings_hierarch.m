function [TCellnew,CommCellnew] = organize_clusterings_hierarch(TCell,CommCell)

% Combine clusterings together
Tlabels = TCell{1};
Comms = CommCell{1};
ssz = size(TCell,2);
for i = 2:ssz
    Tlabels = [Tlabels,TCell{i}];
    Comms = [Comms,CommCell{i}];
end
% sort in increasing order according to the tree node label indices
[Tlabels,I] = sort(Tlabels);
Comms = Comms(I);

TCellnew = cell(1,ssz);
CommCellnew = cell(1,ssz);

for i = 1:(ssz-1)
    [~,CommsTop,TlabelsTop]=nested_filter_int(Comms,Tlabels);
    [CommsTop,TlabelsTop]=eliminate_duplicates_int(CommsTop,TlabelsTop);
    TCellnew{i} = TlabelsTop;
    CommCellnew{i} = CommsTop;
    locs = ismember(Tlabels,TlabelsTop);
    Tlabels(locs) = [];
    Comms(locs) = [];
end

[~,Comms,Tlabels]=nested_filter_int(Comms,Tlabels);
[Comms,Tlabels]=eliminate_duplicates_int(Comms,Tlabels);
TCellnew{ssz} = Tlabels;
CommCellnew{ssz} = Comms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           BELOW IS CALLED BY ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nestnumber,CommunityCell,Tlabels]=nested_filter_int(CommunityCell,Tlabels)
%
sz=size(CommunityCell,2);
%
[~,sizelist]=cellfun(@size,CommunityCell);
%
sizelist=[sizelist',(1:sz)'];
sizelist=sortrows(sizelist);
%
nestnumber=0;
nestpairs=zeros(nchoosek(sz,2),2);
%
for i=1:sz-1;
    bgn=find(sizelist(:,1)>sizelist(i,1),1);
    if isempty(bgn)
        break
    end
    c1=CommunityCell{sizelist(i,2)};
    for j=bgn:sz
        c2=CommunityCell{sizelist(j,2)};
        if (isempty(setdiff(c1,c2))&&~isempty(setdiff(c2,c1)))
            nestnumber=nestnumber+1;
            nestpairs(nestnumber,:)=[sizelist(i,2),sizelist(j,2)];
            %
        end
    end
end
%
nestpairs(nestnumber+1:end,:)=[];
%
%
CommunityCell(nestpairs(:,1)')=[];
Tlabels(:,nestpairs(:,1)')=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CommunityCell,Tlabels]=eliminate_duplicates_int(CommunityCell,Tlabels)

% this program eliminates any duplicate clusters in case they should result
% from any other post processing step

sizelist=cellfun('length',CommunityCell);
% We need to make sure that we do not sort sizelist

vals=unique(sizelist);
szvals=size(vals,2);

for i=1:szvals
    locs=find(sizelist==vals(i));
    szlocs=size(locs,2);
    j=1;
    while j<szlocs
        % get a cluster from the list described by locs
        c1=CommunityCell{locs(j)};
        h=j+1;
        % we can update j now since it will not be needed again until the
        % while loop conditioned on h seen below
        j=j+1;
        while h<=szlocs
            c2=CommunityCell{locs(h)};
            if isequal(c1,c2)
                % here we get rid of the cluster c2 because it is later in
                % the ordering of CommunityCell
                CommunityCell(locs(h))=[];
                Tlabels(locs(h))=[];
                sizelist(locs(h))=[];
                %
                szlocs=szlocs-1;
                locs(h)=[];
                locs(h:end)=locs(h:end)-1;                
            else
                h=h+1;
            end
        end        
    end   
end

