for c=1:length(S.sOCPR)
if S.sOCPR(c,2)== S.eOCPR(c,2)
    startTime = S.sOCPR(c,1);
    endTime = S.eOCPR(c,1);
    timeBetween_p = startTime <= S.period(1:end,1) & S.period(1:end,1) <= endTime;
    timeBetween_d = startTime <= S.decision(1:end,1) & S.decision(1:end,1) <= endTime;
    [tp_r, tp_c] = find(timeBetween_p)
    td_r = find(timeBetween_d)

Str=struct;
%% start
    Time = startTime;
    Trial = S.sOCPR(c,2);
    Period = "1";
    Event = "start";
 

%% dicision
if ~isempty(td_r)    
    Time = unique(S.decision(td_r));
    Trial = Trial(end-1);
    Period = Period(end-1);
    Event = "decision";

    %% period
    Time = S.period(tp_r,1);
    Trial = Trial(end-1);
    Period = extract(S.period(tp_r,2), digitsPattern);
    Event = "start";


    %% end
    Time = endTime;
    Trial = Trial(end-1);
    Period = Period(end-1);
    Event = "end";


%% dicision vs period

P = S.period(tp_r,1);
D = unique(S.decision(td_r));
Time=[];
for i=1:length(P)
    for j=1:length(D)
        if P(i) > D(j)
            Time= [Time; D(j)];
            D(j) = NaN;
        end
    end
    Time = [Time; P(i)];
end
if length(P) < length(D)
    for a=length(P)+1:length(D)
        Time=[Time; D(a)];
    end
elseif length(P) == length(D)
    Time=[Time; D(length(D))];
end





















end
end