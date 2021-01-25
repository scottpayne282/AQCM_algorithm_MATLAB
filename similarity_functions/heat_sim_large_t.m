function [A]=heat_sim_large_t(A,t)

% for G the adjacency data and S the desired nxn similarity matrix
% usage: [S]=heat_sim_large(G);

n=size(A,1);
% A is the adjacency data
A=full(A);
% compute the degree of each node
D=sum(A,2);
% divide each row by degree
for i=1:n
    A(i,:)=A(i,:)/D(i);
end



I=eye(n,n);
A=expm(t*(A-I))-(exp(-t)*I);

S1=pdist(A,'cosine');
S1=squareform(S1);
S1=1-S1;
%
S2=pdist(A','cosine');
S2=squareform(S2);
S2=1-S2;
%
A=(S1+S2)/2;












