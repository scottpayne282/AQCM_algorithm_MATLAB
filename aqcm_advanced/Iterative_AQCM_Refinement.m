function [MainCell]=Iterative_AQCM_Refinement(Comms,Sc,depth,iter)

% Sc should be the contracted graph of S where Comms was the clustering
% used to make Sc
% depth is the number of clusterings we will get in the multi edge cut
% community detection

sz = size(Comms,2);
MainCell = cell(sz,iter+1);
MainCell(:,1) = Comms';

tol=0.008;
lam=1;tee=1;

for jj = 1:iter
    [~,ClustCell,InfoCell,n]=AQCMpar(Sc,tol,lam,tee,5e+25,3);
    %
    [T,TCell,CommCell]=community_detectMINavMemSparseMulti(ClustCell,InfoCell,n,depth);
    clsizes = zeros(1,length(TCell));
    for ii = 1:length(TCell)
        clsizes(ii) = length(TCell{ii});
    end
    choice = find(clsizes==max(clsizes),1);
    Tlabels = TCell{choice};
    Comms = CommCell{choice,:};
    [~,Comms,Tlabels]=nested_filter(Comms,Tlabels);
    %
    [Comms,Tlabels]=post_mpfs(Sc,Comms,Tlabels,T,0);
    [Comms,Tlabels]=eliminate_duplicates(Comms,Tlabels);
    %
    [Comms]=eliminate_multimembers(Sc,Comms,n);
    sizelist = comm_sizelist(Comms)';
    locs = sizelist==0;
    Comms(locs) = [];
    Tlabels(locs) = [];
    [Comms,Tlabels]=eliminate_duplicates(Comms,Tlabels);
    [Tlabels,Comms,~]=addin_singletons(T,Tlabels,Comms,n);
    [Sc]=build_contracted_qcmpPar(Comms,n,Sc,5e+25,3);
    MainCell(1:length(Comms),jj+1) = Comms';
end

[~,ClustCell,InfoCell,n]=AQCMpar(Sc,tol,lam,tee,5e+25,3);
sz2 = size(ClustCell,2);
sz1 = size(ClustCell,1);
NewCell = cell(sz,sz2);
NewCell(1:sz1,:) = ClustCell;
MainCell = [MainCell,NewCell];

