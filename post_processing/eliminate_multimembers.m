function [CommunityCell2]=eliminate_multimembers(S,CommunityCell,n)


sz=size(CommunityCell,2);
C=zeros(n,sz);
for i=1:sz
    labels=CommunityCell{i};
    C(labels,i)=1;
end
intersections=sum(C,2);
%maxoverlap=max(intersections);
numMultimembers=sum(intersections>1);
%
multimembers=find(intersections>1);
%multimemdata=[multimembers,zeros(numMultimembers,maxoverlap)];
Removals=cell(numMultimembers,1);
for i=1:numMultimembers
    v=multimembers(i);
    clusts=find(C(v,:));
    length=size(clusts,2);
    %multimemdata(i,2:(1+length))=clusts;
    %
    contr=zeros(1,length);
    for j=1:length
        labels=CommunityCell{clusts(j)};
        % in the line below we could adjust our calculation by replacing v
        % with the vector multimembers'....that way we would calculate
        % relationships from v to non-multimembers only
        labels=setdiff(labels,multimembers');
        contr(j)=sum(S(v,labels))/(size(labels,2));
    end
    maxcontr=max(contr);
    %using the following code we convert the row array "clusts" of cluster
    %labels to which v belongs into a row array of cluster labels that v
    %wants to be removed from (it no longer wants to be a mamber because
    %its contribution there is not maximum among its various contributions)
    clusts(find(contr==maxcontr))=[];
    Removals{i}=clusts;
end
%
CommunityCell2=CommunityCell;
% the cell "Removals" stores the information about which clusters we want
% to remove each multimember from as a row array of cluster labels, the ith
% cell represents multimember i
for i=1:numMultimembers
   clusts=Removals{i};
   k=size(clusts,2);
   v=multimembers(i);
   for j=1:k
       newcluster=CommunityCell2{clusts(j)};
       newcluster=setdiff(newcluster,v);
       CommunityCell2{clusts(j)}=newcluster;
   end
     
end
%------------------------------------------------------------------------









