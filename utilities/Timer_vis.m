function Timer_vis(Timer)


labelcell={'seeds','growth','merge','contract','update'};
figure;
len=size(Timer,2);
hold on
sz=size(Timer,1);
for i=1:sz
    plot(1:len,Timer(i,:));
end
legend(labelcell);
hold off
f=gcf;
f.Position=[55.4000 77 1396 684.8000];