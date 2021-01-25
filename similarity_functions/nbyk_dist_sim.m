function [S]=nbyk_dist_sim(P)

% the data P should be in the format where a row is a point and the columns
% represent the variables/dimensions

S=pdist(P);
S=squareform(S);
maxS=max(max(S));
S=S-maxS;
S=-S;
maxS=max(max(S));
S=S./maxS;