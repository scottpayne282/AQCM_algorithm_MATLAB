function [clusters]=find_datapoint(CommunityCell,point)


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