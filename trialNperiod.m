
Str=struct('Time',[],'Trial',[],'Period',[],'Event',[]);

var_start = ["scont", "sOCPR"];
var_end = ["endcont", "eOCPR"];

prompt = "OCPR? control? or both? ";
x = input(prompt)
if x == "OCPR"
    var_start = ["sOCPR"];
    var_end = ["eOCPR"];
elseif x == "control"
    var_start = ["scont"];
    var_end = ["endcont"];
elseif x == "both" | isempty(prompt)
    var_start = ["scont", "sOCPR"];
    var_end = ["endcont", "eOCPR"];
else
    error('Input must be one of the three above')
end

for c=1:length(S.scont)
    for e=1:length(var_start)
        A=S.(var_start{e}); B=S.(var_end{e});

    startTime = A(c,1);
    endTime = B(c,1);
    timeBetween_p = startTime <= S.period(1:end,1) & S.period(1:end,1) <= endTime;
    timeBetween_d = startTime <= S.decision(1:end,1) & S.decision(1:end,1) <= endTime;
    [tp_r, tp_c] = find(timeBetween_p);
    td_r = find(timeBetween_d);

%% start
    S_Time = startTime;
    S_Trial = A(c,2);
    S_Period = "1";
    S_Event = "start";

%% make struct

if isempty(Str)
    Str.Time= S_Time;
    Str.Trial= S_Trial;
    Str.Period= S_Period;
    Str.Event= S_Event;
else

    Str.Time= [Str.Time; S_Time];
    Str.Trial= [Str.Trial; S_Trial];
    Str.Period= [Str.Period; S_Period];
    Str.Event= [Str.Event; S_Event];
end
%% dicision
% if ~isempty(td_r)    
    D_Time = unique(S.decision(td_r));
    D_Trial = Str.Trial(end);
    D_Period = Str.Period(end);
    D_Event = "decision";

    %% period
    P_Time = S.period(tp_r,1);
    P_Trial = Str.Trial(end);
    P_Period = extract(S.period(tp_r,2), digitsPattern);
    P_Event = "start";



%% dicision vs period

for i=1:length(P_Time)
    for j=1:length(D_Time)
        if P_Time(i) > D_Time(j)
            Str.Time= [Str.Time; D_Time(j)];
%             D_Time(j) = NaN;
            Str.Trial = [Str.Trial; D_Trial];
            Str.Period = [Str.Period; D_Period];
            Str.Event = [Str.Event; D_Event];
        else 
            P_Time(i) < D_Time(j);
            break;
        end
    end
    Str.Time = [Str.Time; P_Time(i)];
            Str.Trial = [Str.Trial; P_Trial];
            Str.Period = [Str.Period; P_Period(i)];
            Str.Event = [Str.Event; P_Event];
end

if length(P_Time) < length(D_Time)
    for a=length(P_Time)+1:length(D_Time)
        Str.Time=[Str.Time; D_Time(a)];
            Str.Trial = [Str.Trial; D_Trial];
            Str.Period = [Str.Period; D_Period ];
            Str.Event = [Str.Event; D_Event];
    end
elseif length(P_Time) == length(D_Time)
    Str.Time=[Str.Time; D_Time(length(D_Time))];
            Str.Trial = [Str.Trial; D_Trial];
            Str.Period = [Str.Period; D_Period ];
            Str.Event = [Str.Event; D_Event];
end

    %% end
    E_Time = endTime;
    E_Trial = Str.Trial(end);
    E_Period = Str.Period(end);
    E_Event = "end";

    Str.Time= [Str.Time; E_Time];
    Str.Trial= [Str.Trial; E_Trial];
    Str.Period= [Str.Period; E_Period];
    Str.Event= [Str.Event; E_Event];

    end
end

%% make a table
trialNperiod = table(Str.Time, Str.Trial, Str.Period, Str.Event);


