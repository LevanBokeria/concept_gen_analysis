function save_table_for_excel(table_to_write,filename,cleanCellArrays)

% Description: a helper function for writing data tables as either xlsx or
% csv formats. It optionally cleans the columns that have cell arrays,
% since those cannot be written in a csv or excel format.

table_worked_on = table_to_write;

nFields = numel(table_worked_on.Properties.VariableNames);

for iF = 1:nFields
    
    curr_field = table_to_write.Properties.VariableNames{iF};
    
    if cleanCellArrays
        % If this field has cell arrays of length more than 1, delete it
        if isa(table_to_write.(curr_field)(1),'cell')
            if isa(table_to_write.(curr_field){1},'cell')
                
                if numel(table_to_write.(curr_field){1}) > 1
                    
                    table_worked_on.(curr_field) = [];
                end
            end
        end
    end
end % for iF

%% Write the table
writetable(table_worked_on,filename);