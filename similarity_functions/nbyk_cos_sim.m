function [S]=nbyk_cos_sim(P)


oldway = false;

if oldway

    sz=size(P,1);

    S=zeros(sz,sz);
    for i=1:sz
        for j=i:sz
            vi=P(i,:);
            vj=P(j,:);
            lengthvi=sqrt(sum(vi.^2));
            lengthvj=sqrt(sum(vj.^2));
            S(i,j)=(vi*vj')/(lengthvi*lengthvj);
        end
    end
    S1=triu(S,1);
    S=S+S1';
else
    S=pdist(P,'cosine');
    S=squareform(S);
    S=1-S;

end