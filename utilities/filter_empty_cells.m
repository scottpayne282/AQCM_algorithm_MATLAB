function [CommunityCell2,Tlabels2]=filter_empty_cells(CommunityCell2,Tlabels)

% remove cells that represent empty clusters (all were multimembers and all
% left that cluster to join another)
locs=find(cellfun('length',CommunityCell2)==0);

Tlabels2=Tlabels;

Tlabels2(locs)=[];
CommunityCell2(locs)=[];
