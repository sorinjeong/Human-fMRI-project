% Data example
Responses = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1]; % 0 : incorrect trial, 1: correct trial
[pdata cdata] = LearningCurve_WinBugs_230829(Responses);

% Plot example
figure(1);
plotI(Responses, ones(1,length(Responses))); hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'k-'); hold on;
plot(pdata(:,3),'r-'); hold on;
plot(pdata(:,4),'k-');
l=line([cdata cdata],[0 1]); l.Color='b'; l.LineWidth=1;


