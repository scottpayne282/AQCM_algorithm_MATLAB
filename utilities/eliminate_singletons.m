function [CommunityCell2,Tlabels2]=eliminate_singletons(CommunityCell,Tlabels)

% remove cells that represent "size one" clusters 

locs=find(cellfun('length',CommunityCell)==1);

Tlabels2=Tlabels;
CommunityCell2=CommunityCell;

Tlabels2(locs)=[];
CommunityCell2(locs)=[];