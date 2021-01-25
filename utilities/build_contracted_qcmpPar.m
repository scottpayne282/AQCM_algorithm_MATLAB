function [A]=build_contracted_qcmpPar(Comms,n,S,parthresh,wrk)

N=cellfun('length',Comms);

sz=size(Comms,2);
F=false(n,sz);

for i=1:sz
    F(Comms{i},i)=true;
    
end

% sort the communities so that the singletons (communities size 1) are the
% last in the list of communities
[N,sortind]=sort(N,'descend');
F=F(:,sortind);
% sortind is the sort index variable we will use to "unsort" after we build
% the contracted graph A
ind=find(N==1);
grnsetlbls=zeros(1,sz);
for i=1:size(ind,2)
    grnsetlbls(ind(i))=find(any(F(:,ind(i)),2));
end
[~,sortind2]=sort(grnsetlbls);
F=F(:,sortind2);

ind=find(N==1);

if length(ind)<2
    n=size(S,1);
    s=size(F,2);
    numinds=nchoosek(s,2);
    IND = (1:numinds)';
    J = round(floor(-.5 + .5 * sqrt(1 + 8 * (IND - 1))) + 2);
    I = round(J .* (3 - J) / 2 + IND - 1);
    S(1:n+1:n^2)=0;
    
    A=zeros(s,s);
    
    if numinds<parthresh
        for i=1:numinds
            A(I(i),J(i)) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
       -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
       / ...     
    (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  ); 
        end
    else
        %parpool(3);
        k = floor(numinds/wrk);
        m=floor(numinds/k);
        remaind=mod(numinds,k);
        Vals=zeros(k,m);
        parfor ii=1:m
            indexes=(1:k)+k*(ii-1);
            v=zeros(k,1);
            for j=1:k
                i=indexes(j);
                v(j) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
               -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
               / ...     
            (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  );
            end
            Vals(:,ii)=v;
        end
        Vals=[Vals(:);zeros(remaind,1)];
        indexes=(numinds-remaind+1):numinds;
        for j=1:remaind
            i=indexes(j);
            Vals(i) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
           -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
           / ...     
        (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  );
        end
        lininds=sub2ind([s,s],I,J);
        A(lininds)=Vals;
        %delete(gcp('nocreate'))
    end
    A=A+A';
else
    ind=ind(1);
    n=size(S,1);
    s=size(F,2);
    numinds=nchoosek(s,2);
    IND = (1:numinds)';
    J = round(floor(-.5 + .5 * sqrt(1 + 8 * (IND - 1))) + 2);
    I = round(J .* (3 - J) / 2 + IND - 1);
    % get rid of the indices that correspond to singleton clusters
    sinds=find(I>=ind)';
    I(sinds)=[];
    J(sinds)=[];
    numinds=size(I,1);
    
    S(1:n+1:n^2)=0;

    A=zeros(s,s);
    
    if numinds<parthresh
        for i=1:numinds
            A(I(i),J(i)) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
       -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
       / ...     
    (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  ); 
        end
    else
        %parpool(3);
        k = floor(numinds/wrk);
        m=floor(numinds/k);
        remaind=mod(numinds,k);
        Vals=zeros(k,m);
        parfor ii=1:m
            indexes=(1:k)+k*(ii-1);
            v=zeros(k,1);
            for j=1:k
                i=indexes(j);
                v(j) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
               -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
               / ...     
            (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  );
            end
            Vals(:,ii)=v;
        end
        Vals=[Vals(:);zeros(remaind,1)];
        indexes=(numinds-remaind+1):numinds;
        for j=1:remaind
            i=indexes(j);
            Vals(i) =(  sum(sum(  S( F(:,I(i)) , F(:,J(i)) )  )) ...
           -0.5*( sum(sum( S( F(:,I(i)) & F(:,J(i)) , F(:,I(i)) & F(:,J(i)) ) )) )  )...
           / ...     
        (  sum(F(:,I(i)) & F(:,J(i)))*(sum(F(:,I(i)) & F(:,J(i)))-1)/2 + sum(F(:,I(i)) & F(:,J(i)))*(sum( F(:,I(i)) & ~F(:,J(i)) ) + ...
sum( ~F(:,I(i)) & F(:,J(i)) )) + sum( F(:,I(i)) & ~F(:,J(i)) )*sum( ~F(:,I(i)) & F(:,J(i)) )  );
        end
        lininds=sub2ind([s,s],I,J);
        A(lininds)=Vals;
        %delete(gcp('nocreate'))
    end
    %
    list=(find(any(F(:,ind:s),2)))';
    A(ind:s,ind:s)=triu(S(list,list),1);

    A=A+A';
end

% Now unsort the indexing of A
[~,sortind2]=sort(sortind2);
A=A(sortind2,sortind2);
[~,sortind]=sort(sortind);
A=A(sortind,sortind);





