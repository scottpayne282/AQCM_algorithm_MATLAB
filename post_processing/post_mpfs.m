function [CommunityCell,Tlabels]=post_mpfs(S,CommunityCell,Tlabels,T,thresh)

% NOTE about the architecture of this code: This implementation is actually
% a special case of a more general algorithm designed to apply to any
% points in a data set not just unclustered points. For that reason some of
% the code lines will seem unecessary for the case where we work on
% unclustered data only. The code was originally written for the general
% case with runtime speed in mind, that code was then adapted for the use
% seen here but we didn't re-work it completely.


% convert CommunityCell clustering storage to logical array format
n=size(S,1);
sz=size(CommunityCell,2);
F=false(n,sz);
for i=1:sz
    F(CommunityCell{i},i)=true;    
end
%-------------------------------------------------------------------------
% add to clustering any data points not covered yet
Q=diag(~any(F,2));
if any(any(Q,2))
    Q(:,~any(Q,1))=[];
    % working list of data points is the data that was not yet clustered
    list=find(any(Q,2))';
    unclusteredData=true;
    F=[F,Q];  
else
    unclusteredData=false;
end
%------------------------------------------------------------------------

if unclusteredData
    sz=size(F,2);

    k=size(list,2);
    Singlabels=zeros(1,k);
    for i=1:k
        Singlabels(i)=find(T(:,list(i)),1);
    end
    Tlabels=[Tlabels,Singlabels];
    %-------------------------------------------------------------------------
    S(1:n+1:n^2)=0;  % this helps speed up density calculations
    %-------------------------------------------------------------------------
    % calculate cluster densities
    Densities=zeros(1,sz);
    for i=1:sz
        sizec=sum(F(:,i));
        if sizec>1
            Densities(i)=sum(sum(S(F(:,i),F(:,i))))/((sizec)*(sizec-1));
        else
            Densities(i)=1;        
        end
    end
    %-------------------------------------------------------------------------
    %-------------------------------------------------------------------------
    % compute mutual preference cluster choice for each covered node
    k=size(list,2);
    choicelist=zeros(k,1);
    maxpreflist=zeros(k,1);
    for j=1:k
        single=list(j);    
        weightsum=zeros(1,sz);
        weightavg=zeros(1,sz);
        %--------------------------------------------------------------------
        locsj=false(1,sz);
        %--------------------------------------------------------------------
        for i=1:sz
            labels=F(:,i);
            %----------------------------------------------------------------
            if labels(single) % checks if the single datapoint is in cluster i
                locsj(i)=true;
                labels(single)=false;
            end
            %----------------------------------------------------------------
            weightsum(i)=sum(S(single,labels));
            weightavg(i)=weightsum(i)/sum(labels);      
        end
        maxwtav=max(weightavg);
        %---------------------------------------------------------------------
        jDensities=Densities;   
        % here we adjust density values, see explanation below
        locsj=find(locsj);
        for jj=1:size(locsj,2)
            sizec=sum(F(:,locsj(jj)));
            if sizec>2
                jDensities(locsj(jj))=((Densities(locsj(jj))*nchoosek(sizec,2))-...
                    weightsum(locsj(jj)))/nchoosek(sizec-1,2);
            end
        end
        %---------------------------------------------------------------------
        pref=(weightavg.^2)./jDensities; %use this line if we are dividing by 
        % densities of the clusters setminus the data point "single". Otherwise
        % use the code below instead.
        %pref=(weightavg.^2)./Densities;
        pref=pref./maxwtav;
        %
        maxpreflist(j)=max(pref);
        ind=find(pref==maxpreflist(j));
        % the purpose of the following conditional is so that any data point
        % with more than one cluster of maximum preference will have choice
        % "0", this way we may identify data points that are extremely
        % uncertain of which cluster in which they belong
        if size(ind,2)==1
            choicelist(j)=ind;
        end
    end
    %-------------------------------------------------------------------------
    % make a list of datapoints that have a distinct mutual preference....that
    % is, the choicelist value is not zero
    list2=list(choicelist~=0);
    maxpreflist=maxpreflist(choicelist~=0);
    choicelist=choicelist(choicelist~=0);

    %---------------------------------------------------------------------
    % for each point in "list2" (so that choicelist(i)~=0 and its maximum
    % mutual preference value is at least the threshold) join it to its choice
    % cluster
    for i=1:size(list2,2)
        if maxpreflist(i)>thresh
            F(list2(i), choicelist(i))= true;
        end
    end

    %---------------------------------------------------------------------
    % filter out any clusters that became empty
    empties=~any(F,1);
    F(:,empties)=[];
    Tlabels(empties)=[];
    sz=size(F,2);
    % end filter empty clusters
    %---------------------------------------------------------------------
    % filter out any clusters nested in clusters of greater index
    sizelist=[sum(F,1)',(1:sz)'];
    sizelist=sortrows(sizelist);
    nestnumber=0;
    nestpairs=zeros(nchoosek(sz,2),2);
    for i=1:sz-1
        bgn=find(sizelist(:,1)>sizelist(i,1),1);
        if isempty(bgn)
            break
        end
        for j=bgn:sz
            if ~any( F(:,sizelist(i,2)) & ~F(:,sizelist(j,2))  )
                nestnumber=nestnumber+1;
                nestpairs(nestnumber,:)=[sizelist(i,2),sizelist(j,2)];
            end
        end
    end
    nestpairs(nestnumber+1:end,:)=[];
    F(:,nestpairs(:,1)')=[];
    Tlabels(:,nestpairs(:,1)')=[];
    % end filter nested clusters
    %-------------------------------------------------------------------------
    % filter out any singleton clusters
    singletons=find(sum(F,1)==1);
    F(:,singletons)=[];
    Tlabels(singletons)=[];
    % end filter singletons
    %-------------------------------------------------------------------------
    % filter out duplicate clusters (get rid of those with higher index
    sizelist=sum(F,1);
    % We need to make sure that we do not sort sizelist
    vals=unique(sizelist);
    szvals=size(vals,2);
    for i=1:szvals
        locs=find(sizelist==vals(i));
        szlocs=size(locs,2);
        j=1;
        while j<szlocs
            % get a cluster from the list described by locs
            c1=F(:,locs(j));
            h=j+1;
            % we can update j now since it will not be needed again until the
            % while loop conditioned on h seen below
            j=j+1;
            while h<=szlocs
                if isequal(c1,F(:,locs(h)))
                    % here we get rid of the cluster c2 because it is later in
                    % the ordering of CommunityCell
                    F(:,locs(h))=[];
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
    % end filter duplicate clusters
    %---------------------------------------------------------------------
    % make the output variable version of CommunityCell
    sz=size(F,2);
    CommunityCell=cell(1,sz);
    for i=1:sz
        CommunityCell{i}=find(F(:,i))';
    end

end

