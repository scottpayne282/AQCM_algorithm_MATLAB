

%%
% This script demonstrates how to time the individual subroutines and
% visualize their relative runtimes. The functionality is part of the
% AQCMpar version, but we can use it without actually running in parallel
% as demonstrated below.
%
% For this example script we suggest using the larger example data (3824 nodes):
% Load the example data file: network_adjacency_data_3824_nodes.mat
% (it is in the folder: examples)

%% compute the similarity matrix
[S]=heat_sim_large(G);
minS=min(min(S));
S=S-minS;

%% build the hierarchical structure and view runtimes
begintime=clock;
% In the input fields below the value 1e+25 (10^25) is the parallel threshold
% so by choosing such a high value we ensure that parallel processing will
% not kick in (in this example we just want to take advantage of the subroutine
% timer)
[Timer,ClustCell,InfoCell,n]=AQCMpar(S,0.008,1,1,1e+25,3);
% Timer in the output field contains the subroutine runtimes.
endtime=clock;
[qcmtime]=get_qcm_time(begintime,endtime);
% visualize runtimes
Timer_vis(Timer);
% the verticle axis are the runtimes in seconds, the horizontal axis is the
% iteration (level) number...in this example we can see that the
% contraction step takes up most of the runtime overall.





