function [S]=heatsim_NotStrongConnected2(G,t,k,fact)


n=size(G,1);

Gmark=fact*mean(mean(G(G>0)))*ones(n+k,n+k);
L=logical(eye(n+k,n+k));
Gmark(L)=0;
Gmark(1:n,1:n)=G;

[S]=heatsim_w2DiGraph(markovtrans(Gmark),t);
S(:,n+1:n+k)=[];
S(n+1:n+k,:)=[];
