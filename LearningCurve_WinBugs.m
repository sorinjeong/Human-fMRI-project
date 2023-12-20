function [pdata, cback] = LearningCurve_WinBugs(Responses,conf_interval,Chance)
%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2004
%run this script within Matlab
%Modified by Jin Seungwoo 180906

%% Input data
% Responses : behavior data, binary format
% BackgroundProb : chance level, default : 0.5
% conf_interval : confidence interval, default : 95

%% Output data
% pdata : estimate learning curve & confidence interval of upper 95% and lower 5%
% cdata : learning trial

%% Set Path
path_code= 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\code\Learning Analysis';
addpath(genpath([path_code '\Learning Analysis']));
addpath(genpath([path_code '\matbugs-master']));
%%

if ~exist('Responses','var'),  error('Responses is missing!'); end
if ~exist('conf_interval','var'), conf_interval=0.95; end %confience interval
if nargin == 2; BackgroundProb=0.5; disp('Use default: chance level: 0.5');elseif nargin == 3; BackgroundProb = Chance; end % Chance level
fprintf('Use default: conf_interval: %.0f%%\n', conf_interval*100);


MaxResponses=1*ones(1,length(Responses));
dataStruct = struct('n', Responses, 'T', length(Responses), 'ntot', MaxResponses, 'startp', BackgroundProb);

%initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 chains
init1 = struct( 'tau', 0.5);
init2 = struct( 'tau', 1);
init3 = struct( 'tau', 15);

initStructs(1) =  init1;
initStructs(2) =  init2;
initStructs(3) =  init3;

directory = fullfile(path_code,'Learning Analysis','IndividualAnalysisWinBUGS');

[samples, stats, structArray] = matbugs(dataStruct, ...
        fullfile(directory,'Model_bin.txt'), ...
        'init', initStructs, ...
                'nChains', 3, ...
        'view', 0, 'nburnin', 1000, 'nsamples', 5000, ...
        'thin', 10, ...
        'monitorParams', {'p','x','tau','tauprior','sigesq'}, ...
                'Bugdir', [path_code '\WinBUGS14']);

%plot Pr(correct response) (Figure 1 in Smith, Frank, et al. 2004)%%%%%%%%%

pdata =[];
for t = 1:length(Responses) 
        allsamples   = [samples.p(1,:,t) samples.p(2,:,t) samples.p(3,:,t)];
        sort_samples = sort(allsamples);
        total        = length(sort_samples);
        
     if conf_interval == 0.95
        ll           = sort_samples(fix(0.05*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(0.95*total));
     else
        ll           = sort_samples(fix((1-conf_interval)*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(conf_interval*total));         
     end
        
        pdata = [pdata; t ll ml ul];
end

cback = find(pdata(:,2) < BackgroundProb);
if(~isempty(cback))
  if(cback(end) < size(Responses,2) )
       cback = cback(end);
  else
       cback = NaN;
  end
else
  cback = NaN;
end


