function [S]=heatsim_w2DiGraph(W,t)
% this function requires an input W which is the Markov transition matrix
% for the graph/digraph G.
% 
% 
n=size(W,1);
K=expm(t*(W-eye(n,n)))-(exp(-t)*eye(n,n));
%
S1=pdist(K,'cosine');
S1=squareform(S1);
S1=1-S1;
%
S2=pdist(K','cosine');
S2=squareform(S2);
S2=1-S2;
%
S=(S1+S2)/2;