

%%
% This script demonstrates how to run AQCMpar which is the code version
% that parallelizes the contraction step of the algorithm (that is the only
% step that can be parallelized).
% In order to use this code you must have the Matlab Parallel Computing
% Toolbox.
%
% It is important to note that sometimes when the data set is large the
% parallel version of the code can actually take much longer than the
% non-parallel version due to the computing overhead caused by broadcast
% variables. For that reason we consider this version of the code to be
% experimental as we cannot ensure that it is useful for large data.
% 
% For this example script we suggest using the larger example data (3824 nodes):
% Load the example data file: network_adjacency_data_3824_nodes.mat
% (it is in the folder: examples)
%

%% compute the similarity matrix
[S]=heat_sim_large(G);
minS=min(min(S));
S=S-minS;

%% build the hierarchical structure 
% we suggest starting the parallel pool first
poolobj = parpool(3);poolobj.IdleTimeout = 120;
% it is useful to set the IdleTimeout long just in case there are hours
% between AQCM iterations...that way the pool remains active and doesn't
% need to be restarted (which takes time)...we note that long runtimes like
% that have only been observed with large data sets on the order of ~20,000
% data.
%
% NOTE that we use the clock function and out time interpreting function
% instead of using Matlab's tic and toc functions. The reason for this is
% that reportedly Matlab's tic and toc functions do not perform accurately
% when parallel code is running.
begintime=clock;
% In the input fields below the value 100000 is the threshold at which
% parallel computing kicks in...that is, if the number of edges of the
% contracted graph is at least 100000, then it will be done in parallel...
% the user must decide what is the desired threshold.
% The value 3 in the inputs below is the number of parallel
% workers....ideally that value should match the number of workers created
% in the poolobj creation above.
[~,ClustCell,InfoCell,n]=AQCMpar(S,0.008,1,1,100000,3);
% the ~ in the output field is for an output object discussed in the
% example about timing subroutines.
endtime=clock;
[qcmtime]=get_qcm_time(begintime,endtime);


