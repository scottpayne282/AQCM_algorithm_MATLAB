function [X,perm]=dendroLayout(T,n,layout)
% the input "layout" must be the character string 'simple' or the character
% string 'reverse'. The choice 'simple' yields a layout using the
% dendrogram T in a straightforward way. The choice 'reverse' adds a sink
% node to the dendrogram, reverses the direction of each edge and uses this
% new graph to calculate the layered layout

if strcmp(layout,'simple')
    Tree=digraph(T);
    figure('visible','off');
    f = plot(Tree,'Layout','layered');
    %
    X=[(f.XData)',(f.YData)'];
    %
    close gcbf
elseif strcmp(layout,'reverse')
    nodes=size(T,1);
    T=[zeros(1,nodes);T];
    T=[zeros(nodes+1,1),T];
    T(2:n+1,1)=1;
    T=T';
    % now use Matlab's built in layered layout on T
    Tree=digraph(T);
    figure('visible','off');
    f = plot(Tree,'Layout','layered');
    %
    X=[(f.XData)',(f.YData)'];
    %
    close gcbf
    X(:,2)=-X(:,2);
    X(1,:)=[];
else
    error('The input "layout" must be a character string...see function');
end
%---------------------------------
% code to get a permutation from the x data of the leaf row of a dendrogram
x=[X(1:n,1),(1:n)'];
x=sortrows(x);
perm=(x(:,2))';



