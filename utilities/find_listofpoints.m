function [locations]=find_listofpoints(CommunityCell,list)

sz=size(list,1);
locations=zeros(sz,1);
for i=1:sz
    locations(i)=find_datapoint(CommunityCell,list(i));
    
end