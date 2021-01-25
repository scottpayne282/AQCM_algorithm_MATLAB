function [Tsinglabels,CommSingCell,singletons]=addin_singletons(T,Tlabels,...
    CommunityCell,n)
%
% This technical function is designed to get the vertex labels of the data
% points that were not clustered, that is: the points which were not placed
% in any cluster by the community detection algorithm.
% The inputs are data structures that were output by the community
% detection software, n is simply the number of data points of the set we
% are working on.
%------------------------------------------------------------------------
sz=size(CommunityCell,2);
verts=zeros(1,n);
for i=1:sz
    labels=CommunityCell{i};
    verts(labels)=verts(labels)+1;
end
singletons=find(verts==0);
%
sz2=size(singletons,2);
SingCell=cell(1,sz2);
Singlabels=zeros(1,sz2);
for i=1:sz2
    SingCell{i}=singletons(i);
    Singlabels(i)=find(T(:,singletons(i)),1);
end
CommSingCell=[CommunityCell,SingCell];
% now we may update Tlabels
Tsinglabels=[Tlabels,Singlabels];

    