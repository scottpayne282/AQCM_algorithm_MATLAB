

%%
% This script demonstrates how to run AQCM and create a more simplified
% hierarchy tree built from several optimized clusterings selected from the
% AQCM tree. (see the script "aqcm_several_optimized_clusterings.m" for 
% explanations of the multi-clustering code details)
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

%% get k many clusterings
k = 5;
% the value k in the input fields is the number of clusterings desired
[T,TCell,CommCell]=community_detectMINavMemSparseMulti(ClustCell,InfoCell,n,k);
% Now we enforce strict hierarchies (see the script "aqcm_several_optimized_clusterings.m" 
% for explanation)
[TCell,CommCell] = organize_clusterings_hierarch(TCell,CommCell);
% IMPORTANT: TCell{i} is a set of Tlabels of clustering i and CommCell{i}
% is the set of corresponding clusters for clustering i AND the
% organization is that CommCell{1} is the highest hierarchy and the order
% from i=1,2,3... proceeds from higher to lower. THIS ORDERING IS REQUIRED
% IN ORDER FOR THE TREE MAKING CODE BELOW TO WORK.
% 
% We can mention here that the "organize_clusterings_hierarch" function
% also removes any nested clusters so there is no need to run the nested
% cluster removal function that we sometimes use.
%
% We also need to add in any unclustered data at each of the clustering
% levels
for i = 1:k
    [~,Comms,~]=addin_singletons(T,TCell{i},CommCell{i},n);
    CommCell{i} = Comms;
end
%% make a new tree
% make a cellgrid variable representing the k clusterings selected above
[ClustCell] = make_ClustCell_proto(CommCell);
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



