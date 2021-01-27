

%%
% This script demonstrates how to run AQCM and then select a clustering
% "manually" as a specific level of the tree.
% 
% For example data we provide network data as this is tyically what we
% study and we have a useful similarity funtion available.
%
% First load the example data file: network_adjacency_data_1446_nodes.mat
% (it is in the folder: examples)
%
% The example data is a 1446x1446 adjacency matrix from a social media
% network.

%% compute the similarity matrix
[S]=heat_sim_large(G);
minS=min(min(S));
S=S-minS;

%% build the hierarchical structure and convert output to tree format
[ClustCell,InfoCell,n]=AQCM(S,0.008,1,1);
[T,levelsizes,~]=tree_from_cell_grid(ClustCell,n);

%% assign level numbers to tree nodes 
[levels]=assign_tree_levels(T);
% levels is a 1xnn array where nn is the number of nodes of the tree T,
% the ith entry of levels is the level of tree node i
%% find the nodes at level 2 and output clustering
Tlabels = find(levels==2);
[Comms]=leaves_from_tree_nodes(T,n,Tlabels);
% Comms is the clustering obtained from T at level 2



