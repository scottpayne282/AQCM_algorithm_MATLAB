function [Comms]=leaves_from_tree_nodes(T,n,roots)

% assume indexes 1:n is the bottom row (bottom level) of the rooted
% hierarchy tree T (all arrows point down)

% the input variable "roots" is a row array of numeric labels of nodes of
% the tree T for which we wish to find the leaves

% the stack variable stk stores items by row where the first entry is the
% node number of a node and the second is it's level inthe hierarchy tree

num=size(T,1);
sz = size(roots,2);
Comms = cell(1,sz);

for ii=1:sz
    % initialize storage for the leaves
    leaves = false(1,n);
    % initialize a stack
    stk=zeros(num,1);
    hgt=1; % current stack height
    % initialize the stk with the root as the only item
    stk(1)=roots(ii);
    % initialize a variable to keep track of previous children nodes
    prev_chld = [];
    while hgt>0
        % get a node from the top of the stack
        node=stk(hgt);
        stk(hgt)=0;
        hgt=hgt-1;
        % get the children of that node
        chld=find(T(node(1),:));
        nchld=size(chld,2);
        if nchld>0
            % record the current children
            prev_chld = chld;
            % push the children to the stack
            stk((hgt+1):(hgt+nchld))=chld';
            % update the stack height
            hgt=hgt+nchld;
        else
            leaves(prev_chld) = true;
        end
    end
    Comms{ii} = find(leaves);
end
% 





