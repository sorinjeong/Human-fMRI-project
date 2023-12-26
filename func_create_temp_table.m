function temp_table = func_create_temp_table(DMTS, event_table, varTypes)
    temp_table = array2table(cell(height(DMTS), width(event_table)), 'VariableNames', event_table.Properties.VariableNames);

    % 각 변수의 형식을 event_table의 형식과 동일하게 설정합니다.
    for i = 1:width(event_table)
        if strcmp(varTypes{1,i}, 'double')
            temp_table.(i) = NaN(height(DMTS), 1);
        elseif strcmp(varTypes{1,i}, 'cell')
            temp_table.(i) = cell(height(DMTS), 1);
        elseif strcmp(varTypes{1,i}, 'string')
            temp_table.(i) = strings(height(DMTS), 1);
        else
                    % 다른 형식의 경우에 대한 처리를 여기에 추가할 수 있습니다.
        end
    end
end
