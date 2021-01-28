

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

%% get three clusterings
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
