function [lssNames, lssOnsets, lssDurations] = lssMakeVectors(originalNames, originalOnsets, originalDurations, includeConditions)
% FORMAT [lssNames, lssOnsets, lssDurations] = lssMakeVectors(originalNames, originalOnsets, originalDurations, includeConditions)
% Uses SPM-format vectors (in variables corresponding to names, onsets, and
% durations) to create cell arrays of names, onsets, and durations for each
% LS-S model for conditions of interest.
%
%
% input
% originalNames:        Cell array of condition names for a single block.
% originalOnsets:       Cell array of trial onset vectors for a single block.
% originalDurations:    Cell array of trial duration vectors for a single block.
% includeConditions:    A cell array of events we wish to model with the
%                       beta LSS method.
% 
% output
% lssNames:             Cell array of LS-S condition names for a single
%                       block. Format lssNames{includedCondition}{conditionName}
% lssOnsets:            Cell array of LS-S condition onsets for a single
%                       block. Format lssOnsets{includedCondition}{trialVersion}(trial)
% lssDurations:         Cell array of LS-S condition durations for a single
%                       block. Format lssDurations{includedCondition}{trialVersion}(trial)
for iCond = 1:length(includeConditions)
    % Determine where conditions of interest are in vectors.
    % Setdiff reorders conditions, otherConditions must be reordered.
    otherConditionsIdx = ~strcmp(includeConditions{iCond}, originalNames);
    [otherConditions, originalOrder] = setdiff(originalNames, includeConditions{iCond});
    [~, sortedOrder] = sort(originalOrder);
    otherConditions = otherConditions(sortedOrder);
    includeConditionIdx = find(~otherConditionsIdx);
    
    % Check that condition of interest has more than one trial.
    % If condition A only has one trial, you don't need both ConditionA_001
    % and Other_ConditionA, because Other_ConditionA would be empty.
    if ~isempty(setdiff(originalOnsets{includeConditionIdx}, originalOnsets{includeConditionIdx}(1)))
        for jOnset = 1:length(originalOnsets{includeConditionIdx}),
            % Create list of condition names
            % (e.g. ConditionA_001, Other_ConditionA, ConditionB, ConditionC, etc.)
            lssNames{iCond}{jOnset} = [{[originalNames{includeConditionIdx} '_' sprintf('%03d', jOnset)]...
                                        ['Other_' originalNames{includeConditionIdx}]}...
                                       otherConditions];
            
            % Single trial
            lssOnsets{iCond}{jOnset}{1} = originalOnsets{includeConditionIdx}(jOnset);
            lssDurations{iCond}{jOnset}{1} = originalDurations{includeConditionIdx}(jOnset);
            
            % Other trials of same condition
            lssOnsets{iCond}{jOnset}{2} = originalOnsets{includeConditionIdx};
            lssOnsets{iCond}{jOnset}{2}(jOnset) = [];
            lssDurations{iCond}{jOnset}{2} = originalDurations{includeConditionIdx};
            lssDurations{iCond}{jOnset}{2}(jOnset) = [];

            % Other conditions
            counter = 3; % A counter adjusts around the skipped condition.
            for kCond = find(otherConditionsIdx)
                lssOnsets{iCond}{jOnset}{counter} = originalOnsets{kCond};
                lssDurations{iCond}{jOnset}{counter} = originalDurations{kCond};
                counter = counter + 1;
            end
        end
    else
        % Single trial
        lssNames{iCond}{1} = [{[originalNames{includeConditionIdx} '_' sprintf('%03d', 1)]} otherConditions];
        lssOnsets{iCond}{1}{1} = originalOnsets{includeConditionIdx}(1);
        lssDurations{iCond}{1}{1} = originalDurations{includeConditionIdx}(1);
        
        % Other conditions
        conditionCounter = 2; % A counter adjusts around the skipped condition.
        for kCond = find(otherConditionsIdx)
            lssOnsets{iCond}{1}{conditionCounter} = originalOnsets{kCond};
            lssDurations{iCond}{1}{conditionCounter} = originalDurations{kCond};
            conditionCounter = conditionCounter + 1;
        end
    end
end
end