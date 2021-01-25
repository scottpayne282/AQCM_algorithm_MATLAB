function [CommunityCell,Tlabels]=eliminate_duplicates(CommunityCell,Tlabels)

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
