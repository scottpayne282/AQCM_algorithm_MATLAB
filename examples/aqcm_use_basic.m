

%%
% This script demonstrates basic use of AQCM with edge cut cluster
% selection.
% 
% For example data we provide network data as this is tyically what we
% study and we have a useful similarity funtion available.
%
% First load the example data file: network_adjacency_data_3824_nodes.mat
% (it is in the folder: examples)
%
% The example data is a 3824x3824 adjacency matrix from a social media
% network.

%% compute the similarity matrix
[S]=heat_sim_large(G);
% the following is not required but we often do it to zero the minimum
minS=min(min(S));
S=S-minS;

%% build the hierarchical structure
% note that there are some constants in the input fields...we always use
% these
[ClustCell,InfoCell,n]=AQCM(S,0.008,1,1);
% ClustCell describes the structure of the tree:
% ClustCell is a sxk cell array where s is the number of clusters +
% unclustered points at the first level and k is the number of levels.
% ClustCell{i,1} is the numeric labels of the data in cluster i at level 1.
% ClustCell{i,2} is the numeric labels of the level 1 clusters which
% combined to form level 2 cluster i.
% InfoCell contains the information about cluster size and density necessary 
% to compute edge weights for cluster detection. 
% InfoCell{i,j} has the info for the cluster represented by ClustCell{i,j}.
% InfoCell{i,j} is a 2x1 array where the first entry is similarity density
% of the cluster and the second is the number of data represented by the
% cluster (i.e. its size)
% NOTE that in the code below we convert the hierarchy tree info from the
% output structures described above into a standard format in the variable
% T. If only the unweighted tree is desired run: 
%       [T,levelsizes,levels]=tree_from_cell_grid(ClustCell,n);
% Also, n is just the number of data of the underlying data set.

%% weight the tree and select clustering
% now create the weighted tree T as a sparse array and also compute the
% optimal clustering using average edge cut method.
[T,Tlabels,Comms]=community_detectMINavMemSparse(ClustCell,InfoCell,n);
% T(i,j) is the weight of directed edge ij, Tlabels is a row array of node 
% labels of T which represented the selected clustering. Comms is a row cell
% variable such that Comms{i} is a row array of numeric labels of the data
% points of cluster i.
%
% remove potential nested clusters....this can occure somtimes since the
% tree T supports multimembership so edge cuts can possible produce this
% effect though it tends to be not more than a few.
[nestnumber,Comms,Tlabels]=nested_filter(Comms,Tlabels);

%% draw diagrams os the hierarchy tree, clustering, and similarity data
[XT,permT]=dendroLayout(T,n,'simple');
% XT is a mx2 array where row i are the x and y coords of node i of the
% tree T, permT is a permutation of the underlying data set that arranges
% the data according to the layout of the tree T (bottom level of T are the
% data)
%
% DIAGRAM of the tree
dendroOfXColor(T,Tlabels,n,8,XT);
% If a dendrogram with all nodes numerically labeled is desired use:
%           dendroOfX2(T,Tlabels,n,8,XT);
%
% DIAGRAM of the similarity data arranged by the tree
% the number 255 is a resolution for the number of shades for the weighted
% similarity data, 127 or 63 give lower res...1 gives just black or white
% and is a useful setting if the input S is an adjacency matrix
dendro_adjacencyColor(S,permT,255,1,1,T,XT,Tlabels);
%
% NOTE if colored labels on the tree diagrams are not desireable just input
% [] to the functions instead of Tlabels










