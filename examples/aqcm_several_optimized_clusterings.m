

%%
% This script demonstrates how to run AQCM and then get several different
% clusterings by taking successive disjoint edge-cuts, each optimal when
% excluding the previous cuts.
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

%% get two clusterings
% the value 2 in the input fields is the number of clusterings desired
[T,TCell,CommCell]=community_detectMINavMemSparseMulti(ClustCell,InfoCell,n,2);
% TCell contains the Tlabel variables and CommCell contains the
% clusterings, that is CommCell{i} is itself a cell variable that is a
% clustering (see basic example for our standard clustering output vars)

%% compare the two clusterings
% make a layout for the tree T
[XT,permT]=dendroLayout(T,n,'simple');
% show each of the clusterings
Tlabels1 = TCell{1};
dendroOfXColor(T,Tlabels1,n,8,XT);
Tlabels2 = TCell{2};
dendroOfXColor(T,Tlabels2,n,8,XT);

%% enforce hierarchy to clusterings
% as we can see in the tree diagrams created above, neither of the two
% clusterings properly contains the other....that is, the first clustering
% is mostly contained by the second, but there are some clusters in the
% second that are subclusters of clusters in the first.
% We can reorganize the clusters from the two clusterings into two properly
% hierarchical clusterings:
[TCellnew,CommCellnew] = organize_clusterings_hierarch(TCell,CommCell);
% We can mention here that the "organize_clusterings_hierarch" function
% also removes any nested clusters so there is no need to run the nested
% cluster removal function that we sometimes use.
Tlabels1 = TCellnew{1};
dendroOfXColor(T,Tlabels1,n,8,XT);
Tlabels2 = TCellnew{2};
dendroOfXColor(T,Tlabels2,n,8,XT);
% in the above outputs the first set of clusters is strictly higher in the
% hierarchy tree than the second. 


