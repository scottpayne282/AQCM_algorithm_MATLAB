function [sizelist]=comm_sizelist(CommunityCell)

sizelist=cellfun('length',CommunityCell);
sizelist=sizelist';