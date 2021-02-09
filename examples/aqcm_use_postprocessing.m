

%%
%
% This script demonstrates how to run AQCM and then perform simple "postprocessing"
% tasks which can be helpful in certain applications. The post processing
% algorithms can force unclustered data to cluster and also remove
% multi-membership. These two functions are sometimes desirable in data
% analytic jobs.
%
% First load the example data file: network_adjacency_data_1446_nodes.mat
% (it is in the folder: examples)
%

%% compute the similarity matrix
[S]=heat_sim_large(G);
minS=min(min(S));
S=S-minS;

%% build the hierarchical structure
[ClustCell,InfoCell,n]=AQCM(S,0.008,1,1);

%% choose clustering and refine
% here we get the first optimal clustering
[T,Tlabels,Comms]=community_detectMINavMemSparse(ClustCell,InfoCell,n);
[nestnumber,Comms,Tlabels]=nested_filter(Comms,Tlabels);
% The following post processing algorithm causes unclustered data points to
% cluster...depending on the nature of the data the unclustered points may
% prefer to join an existing cluster or they may form clusters with other
% unclustered points....this type of post processing is useful when a
% clustering covering the data set fully is desireable (which is the case here)
% The parameter 0 in the input fields below means that all unclustered
% points are eligible to join clusters (choosing a higher value between zero 
% and 1 means that there is a criteria for a point to have a certain value of
% cluster preference in order to join)
[Comms,Tlabels]=post_mpfs(S,Comms,Tlabels,T,0);
% The following further refines the clustering by eliminating
% multimembership, here this feature is desireable for the purpose of our
% example however it is not a mathematical requirement.
[Comms]=eliminate_multimembers(S,Comms,n);
% the following "clean up" codes take care of the case when some of the
% storage cells in the variable Comms become empty after running the above
% combination of post processing...also we add in any unclustered data that
% are still unaccounted for.
sizelist = comm_sizelist(Comms)';
locs = sizelist==0;
Comms(locs) = [];
Tlabels(locs) = [];
[Tlabels,Comms,~]=addin_singletons(T,Tlabels,Comms,n);


