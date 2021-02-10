function [Timer,ClustCell,InfoCell,n]=AQCMpar(S,tol,lam,tee,parthresh,wrk)
% INPUTS: S is an nxn similarity matrix
%         
%         tol is a numerical tolerance (error) parameter. Experimentally we
%         have found that we prefer to use tol = 0.008, we treat this as a
%         constant although one might experiment by changing value.
%
%         lam, tee were parameters in previous published works but have
%         since been established as lam = 1 and tee = 1
%
%         parthresh is the threshold of the contractionstep job size that
%         will trigger the parallel code to run.....if we want to use this
%         software version without triggering parallel just set parthresh
%         very high like 10^25. We have found that the overhead caused by
%         parallelization sometimes results in a longer runtime than if we
%         just run without parallel.
%
%         wrk is the number of workers desired for parallel....it is best
%         to choose the number of cores you wish to use....for example on a
%         4 core machine we use 3 workers (leaves a core available for the machine)
%
%   This version is the same as AQCM except that we have 
%   parallelized the function contractionstep and added a timer for each of
%   the subroutines
%
%------------------------------------------------------------------------
Timer=zeros(5,200);
% run through first pass
tic;
[Es]=seedselection(S,0.0000000001);
Timer(1,1)=toc;
tic;
[F,D,N]=growthstep(S,Es,tol,lam,tee);
Timer(2,1)=toc;
tic;
[F,D,N]=adjuststep(F,D,N,S,0.0000000001);
Timer(3,1)=toc;

[F,D,N]=addmissedpoints(F,D,N);
%-------------------------------------------------------------------------
% initialize storage structures for the outputs of passes
% pre-allocation is done to improve performance of the storage method 
% when looping through additional passes
n=size(S,1);
s=size(F,2);
ClustCell=cell(s,floor(n/3));
InfoCell=cell(s,floor(n/3));
% store the output of the first pass
for i=1:s;
    ClustCell{i,1}=int16(find(F(:,i)'));
    InfoCell{i,1}=[D(i);N(i)];
end
%-------------------------------------------------------------------------
begintime=clock;
[A]=contractionstep(F,N,S,parthresh,wrk);
endtime=clock;
Timer(4,1)= get_par_time(begintime,endtime);
disp('contraction 1 finished')
%-------------------------------------------------------------------------
% enter the iterative phase of the seed,grow,merge,contract cycle
level=1;
while s>1
    level=level+1;
    %---------------------------------------------------------------------
    tic;
    [Es]=seedselection(A,0.000000000001);
    Timer(1,level)=toc;
    tic;
    [P]=growthstepIterPhase(A,Es,tol,lam,tee);
    Timer(2,level)=toc;
    % any unclustered points of the contracted matrix A will be added into
    % the current clustering F by the update function
    tic;
    [F,P,D,N]=updatestep(F,P,S);
    Timer(5,level)=toc;
    tic;
    [F,P,D,N]=adjuststepIterPhase(F,P,D,N,S,0.0000000001);
    Timer(3,level)=toc;
    %---------------------------------------------------------------------
    snew=size(F,2);
    if snew<s
        s=snew;
    else
        error('number of clusters is at least that of previous iteration');
    end
    % store the output of the seed,grow,merge cycle
    for i=1:s;
        ClustCell{i,level}=int16(find(P(:,i)'));
        InfoCell{i,level}=[D(i);N(i)];
    end
    %---------------------------------------------------------------------
    if s>1
        begintime=clock;
        [A]=contractionstep(F,N,S,parthresh,wrk);
        endtime=clock;
        Timer(4,level)= get_par_time(begintime,endtime);
        disp('contraction finished')
    end
end
% trim off the unused cell blocks that had been pre-allocated but not
% needed
ClustCell(:,level+1:end)=[];
InfoCell(:,level+1:end)=[];
Timer(:,level+1:end)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           BELOW ARE CALLED BY THE ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time]=get_par_time(begintime,endtime)

% this code doesn't give correct if new year rolls over (lol)
% also, doesn't work in February rollover during leap year
% also it assumes month rollover is just to next month, not two months

% the output time is in seconds

% the following is coded to have easy to understand variables in case I
% forget how this works later.
hr1 = begintime(4);
mn1 = begintime(5);
sc1 = begintime(6);
hr2 = endtime(4);
mn2 = endtime(5);
sc2 = endtime(6);

% convert to seconds "past midnight"
t1 = (hr1*60*60) + (mn1*60) + sc1;
t2 = (hr2*60*60) + (mn2*60) + sc2;

% check if month rolled over to next month
mo1 = begintime(2);
mo2 = endtime(2);
dy1 = begintime(3);
dy2 = endtime(3);
if mo2 ~= mo1
    months = [31 28 31 30 31 30 31 31 30 31 30 31];
    dy2 = dy2 + months(mo1);
end

% check if day rolled over to next day
if dy2 ~= dy1
    t2 = t2 + 24*60*60*(dy2-dy1);
end

% just subtract
time = t2 - t1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,D,N]=addmissedpoints(F,D,N)

Q=diag(~any(F,2));
if any(any(Q,2));
    Q(:,~any(Q,1))=[];
    F=[F,Q];
    q=size(Q,2);
    D=[D,Inf(1,q)];
    N=[N,ones(1,q)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A]=contractionstep(F,N,S,parthresh,wrk)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,D,N]=growthstep(S,Es,tol,lam,tee)

n=size(S,1);

F=false(n,n);
D=zeros(1,n);
N=zeros(1,n);
c=false(n,1);
count=1;

while ~isempty(Es)
    %---------------------------------------------------------------------
    % starting a cluster from a seed
    c(Es(1,2:3))=true;
    lininds=sub2ind([n,n],Es(1,2),Es(1,3))';
    Es(1,:)=[];
    sizec=2;
    finishedgrowing=false;
    
    while ~finishedgrowing
        %---------------------------------------------------------------------
        % calculating contributions
        cont=[(1/sizec)*sum(S(:,c),2),(1:n)'];
        cont(c,:)=[];
        maxcont=max(cont(:,1));
        alph=1-(1/(2*lam*(sizec+tee)));
        denC=2*sum(S(lininds))/(sizec*(sizec-1)); 
        
        %---------------------------------------------------------------------
        if maxcont>=alph*denC
            joinups=(cont(cont(:,1)<=maxcont+tol & cont(:,1)>=maxcont-tol,2))'; 
            c(joinups)=true;
            sizec=sizec+size(joinups,2); % update current value cluster size
            temp=nchoosek(find(c),2);
            lininds=sub2ind([n,n],temp(:,1),temp(:,2))';
            locs=ismember(Es(:,4) , lininds); 
            Es(locs,:)=[];
        else
            F(:,count)=c;
            D(count)=2*sum(S(lininds))/(sizec*(sizec-1)); 
            N(count)=sizec;
            count=count+1;
            c(c)=false;
            finishedgrowing=true;
        end
        %----------------------------------------------------------------------
    end
end
%----------------------------------------------------------------------------
F(:,count:end)=[];
D(count:end)=[];
N(count:end)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P]=growthstepIterPhase(A,Es,tol,lam,tee)

n=size(A,1);

P=false(n,n);
c=false(n,1);
count=1;

while ~isempty(Es)
    %---------------------------------------------------------------------
    % starting a cluster from a seed
    c(Es(1,2:3))=true;
    lininds=sub2ind([n,n],Es(1,2),Es(1,3))'; 
    Es(1,:)=[];
    sizec=2;
    finishedgrowing=false;
    
    while ~finishedgrowing
        %---------------------------------------------------------------------
        % calculating contributions
        cont=[(1/sizec)*sum(A(:,c),2),(1:n)'];
        cont(c,:)=[];
        maxcont=max(cont(:,1));
        alph=1-(1/(2*lam*(sizec+tee)));
        denC=2*sum(A(lininds))/(sizec*(sizec-1));
        
        %---------------------------------------------------------------------
        if maxcont>=alph*denC;
            joinups=(cont(cont(:,1)<=maxcont+tol & cont(:,1)>=maxcont-tol,2))'; 
            c(joinups)=true;
            sizec=sizec+size(joinups,2); % update current value cluster size
            temp=nchoosek(find(c),2);
            lininds=sub2ind([n,n],temp(:,1),temp(:,2))';
            locs=ismember(Es(:,4) , lininds);
            Es(locs,:)=[];
        else
            P(:,count)=c;
            count=count+1;
            c(c)=false;
            finishedgrowing=true;
        end
        %----------------------------------------------------------------------
    end
end
%----------------------------------------------------------------------------
P(:,count:end)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,D,N]=adjuststep(F,D,N,S,tol)
%
% sort the clusters in descending order by their density
[D,I]=sort(D,'descend');
F=F(:,I);
N=N(I);
s=size(D,2); % this is the number of clusters
n=size(F,1); % this is the number of data points
S(1:n+1:n*n)=0; % this helps the density calc
j=1;

while j<s
    mergesj=j+find( sum(repmat(F(:,j),[1 s-j]) & F(:,j+1:s)) >...
        0.5*min(repmat(N(j),[1 s-j]),N(j+1:s)) );
    if isempty(mergesj)
        j=j+1;
    else
        m=length(mergesj);
        densities=zeros(1,m);
        for i=1:m;
            % calculate density of the union of clusters j and mergesj(i)
            union=F(:,j) | F(:,mergesj(i));
            sizec=sum(union);
            densities(i)=sum(sum(S(union,union)))/(sizec*(sizec-1));            
        end
        mind=min(densities);
        choice=mergesj(find(densities<=mind+tol,1,'last'));
        
        F(:,choice)=[];
        
        D(choice)=[];
        
        N(choice)=[];
        s=s-1;
        
        %----------------------------------------------------------------
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,P,D,N]=adjuststepIterPhase(F,P,D,N,S,tol)
%
% the following code segment is necessary to keep "single point" clusters
% sorted at the "low end" of the list, this prevents an error that occurs
% in contraction otherwise
D(N==1)=-1;
% sort the clusters in descending order by their density
[D,I]=sort(D,'descend');
F=F(:,I);
P=P(:,I);
N=N(I);
s=size(D,2); % this is the number of clusters
n=size(F,1); % this is the number of data points
S(1:n+1:n*n)=0; % this helps the density calc
j=1;

while j<s
    mergesj=j+find( sum(repmat(F(:,j),[1 s-j]) & F(:,j+1:s)) >...
        0.5*min(repmat(N(j),[1 s-j]),N(j+1:s)) );
    if isempty(mergesj)
        j=j+1;
    else
        m=length(mergesj);
        densities=zeros(1,m);
        for i=1:m;
            % calculate density of the union of clusters j and mergesj(i)
            union=F(:,j) | F(:,mergesj(i));
            sizec=sum(union);            
            densities(i)=sum(sum(S(union,union)))/(sizec*(sizec-1));            
        end
        maxd=max(densities);
        choice=mergesj(find(densities>=maxd-tol,1));
        union=F(:,j) | F(:,choice);
        unionP=P(:,j) | P(:,choice);
        sizec=sum(union);        
        denUnion=sum(sum(S(union,union)))/(sizec*(sizec-1));        
        %----------------------------------------------------------------
        F(:,j)=union;
        F(:,choice)=[];
        P(:,j)=unionP;
        P(:,choice)=[];
        D(j)=denUnion;
        D(choice)=[];
        N(j)=sizec;
        N(choice)=[];
        s=s-1;
        [D,I]=sort(D,'descend');
        F=F(:,I);
        P=P(:,I);
        N=N(I);
        %----------------------------------------------------------------
    end
    
end
% the following code segment resets the density values of single point
% clusters to Inf
D(N==1)=Inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Es]=seedselection(S,tol)

% "tol" is the tolerance for computational error when deciding whether two
% values are equal in the while loop.
% It is recommended to try tol=0.000000000001 as a standard before
% experimenting with other values

n=size(S,1);

if n>2
    S(1:n+1:n*n) = -1;
    [B,I]=sort(S,'descend');
    I=int16(I);
    B(n,:)=[];
    I(n,:)=[];
    Drops=-diff(B);
    Meds=median(Drops,1);

    L=false(n,n);

    for i=1:n
        k=find(Drops(:,i)>=Meds(i),1);
        while B(k,i)-B(k+1,i)<tol;
            k=k+1;
            if k==n-1;
                break
            end
        end
        L(I(1:k,i),i)=true;
    end

    L=triu(L&L',1);
    indx=find(L);
    [row,col]=ind2sub([n,n],indx);
    Es=[S(L),row,col,indx];
    Es=flipud(sortrows(Es));
else
    Es=[S(1,2),1,2,3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,P,D,N]=updatestep(F,P,S)

%-------------------------------------------------------------------------
% add in to P any data points of the contracted matrix A
Q=diag(~any(P,2));
if any(any(Q,2));
    Q(:,~any(Q,1))=[];
    P=[P,Q];
end
%-------------------------------------------------------------------------
n=size(F,1);
S(1:n+1:n*n)=0; % this helps the density calc
s=size(P,2);
F2=false(n,s);
D=zeros(1,s);
N=zeros(1,s);
for i=1:s
    F2(:,i)=any(F(:,P(:,i)'),2);
    N(i)=sum(F2(:,i));
    if N(i)>1;
        D(i)=sum(sum(S(F2(:,i),F2(:,i))))/((N(i))*(N(i)-1));
    else
        D(i)=Inf;
    end
end
F=F2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            END FUNCTIONS

