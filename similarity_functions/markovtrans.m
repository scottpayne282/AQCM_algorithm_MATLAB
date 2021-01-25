function [W]=markovtrans(G)
%
% produce the random walk transition operator matrix W
%the following computes the diagonal "inverse degree matrix"
D=sum(G,2);
D=1./D;
D=diag(D);
%the line below computes the "random walk matrix" for the weighted graph G
W=D*G;