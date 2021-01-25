function [nestnumber,CommunityCell,Tlabels]=nested_filter(CommunityCell,Tlabels)
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

