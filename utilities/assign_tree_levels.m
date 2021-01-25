function [levels]=assign_tree_levels(T)

% assume that the root has the highest index and 1:n is the bottom row

% the stack variable stk stores items by row where the first entry is the
% node number of a node and the second is it's level inthe hierarchy tree

num=size(T,1);
levels=zeros(1,num);
stk=zeros(num,2);
hgt=1; % current stack height
% initialize the stk with the root as the only item
stk(1,1)=num;
stk(1,2)=0; % not necessary, just a reminder that the root level is zero

while hgt>0
    % get a node from the top of the stack
    node=stk(hgt,:);
    stk(hgt,:)=[0 0];
    hgt=hgt-1;
    % get the children of that node
    chld=find(T(node(1),:));
    nchld=size(chld,2);
    if nchld>0
        % push the children to the stack
        stk((hgt+1):(hgt+nchld),1)=chld';
        % store the childrens' levels to the stack
        stk((hgt+1):(hgt+nchld),2)=node(2)-1;
        % store the childrens' levels to the storage variable
        levels(chld)=node(2)-1;
        % update the stack height
        hgt=hgt+nchld;
    end
end

% 
levels=levels-min(levels);




