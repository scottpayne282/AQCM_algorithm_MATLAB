

%%
% IMPORTANT note about this script: This example is a demonstration of how
% to build a pipeline in which AQCM and cluster detection serve as
% subroutines in a larger process. We demonstrate the use of the function 
% Iterative_AQCM_Refinement.m which is itself an experimental
% implementation and is not intended to be a provably optimal approach,
% rather it is one way that one might implement the idea of further iteration.
%
% This script demonstrates how to run AQCM and then create a "contracted graph"
% which will become the input for a further iterative process in which AQCM
% and cluster detection are run subsequently in a similar manner to the
% AQCM tree building process itself. The result of this further iterative
% process is a simplified hierarchy tree built from the subsequent outputs
% of cluster detection...To explain it another way, we are using an
% aglomerative approach to building a hierarchy tree and in our approach
% "AQCM with cluster detection" serves as the subroutine which finds the
% clustering at each level.
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

%% build the contracted graph
% in the following code the input fields 5e+25 and 3 are parameters related
% to using a parallel version, by using the large threshold 5e+25 the
% parallelization does not kick in (see the script aqcm_use_parallel_version.m
% for technical details)
[Sc]=build_contracted_qcmpPar(Comms,n,S,5e+25,3);

%% Input the contracted graph and clustering into repeated iterative process
% In the input fields the parameter 2 controls the number of optimal cuts
% we will try at each level (we end up choosing the one that has the most clusters)
% and the parameter 2 is the number of rounds to iterate.
% WARNING: choosing the iteration number parameter too can cause an error,
% for example with the data set demonstrated here (n=1446) an error happens
% if we ask for 3 or more iterations, but for 2 iterations we get nice
% output.
[ClustCell]=Iterative_AQCM_Refinement(Comms,Sc,3,2);
% the output ClustCell describes a tree (it is the same type of variable output by AQCM)

%% draw the simplified hierarchy tree
% the following "tree-deepener" is for visual enhancement of the tree only 
% and does not affect the nature of the clusterings described
[ClustCell]=tree_deepener(ClustCell,1);
[T,~,~]=tree_from_cell_grid(ClustCell,n); %<-- n is the number of data
% create a layout for T and also apply visual enhancement with "tree_deepener_layout"
[XT,permT]=dendroLayout(T,n,'simple');
[XT]=tree_deepener_layout(XT,0.6);
% draw the diagram (input [] instead of Tlabels since we aren't choosing a 
% specific clustering here)
dendroOfXColor(T,[],n,8,XT);
dendro_adjacencyColor(S,permT,255,1,1,T,XT,[]);
