function [colors]=get_dist_colors(n)


if n>=3
    colors = distinguishable_colors(n+2);
    colors(4:5,:)=[];
else
    colors = distinguishable_colors(n);
    
end