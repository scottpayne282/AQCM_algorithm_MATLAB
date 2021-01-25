function [time]=get_qcm_time(begintime,endtime)

% this code doesn't give correct if new year rolls over (lol)
% also, doesn't work in February rollover during leap year
% also it assumes month rollover is just to next month, not two months

% the output time is in seconds

% the following is coded to have easy to understand variables in case I
% forget how this works later.
hr1 = begintime(4);
mn1 = begintime(5);
sc1 = begintime(6);
hr2 = endtime(4);
mn2 = endtime(5);
sc2 = endtime(6);

% convert to seconds "past midnight"
t1 = (hr1*60*60) + (mn1*60) + sc1;
t2 = (hr2*60*60) + (mn2*60) + sc2;

% check if month rolled over to next month
mo1 = begintime(2);
mo2 = endtime(2);
dy1 = begintime(3);
dy2 = endtime(3);
if mo2 ~= mo1
    months = [31 28 31 30 31 30 31 31 30 31 30 31];
    dy2 = dy2 + months(mo1);
end

% check if day rolled over to next day
if dy2 ~= dy1
    t2 = t2 + 24*60*60*(dy2-dy1);
end

% just subtract
time = t2 - t1;