%% Use this file to make changes to CellList and to the CellParameters Structure of a cell you've already analyzed.
%(along with the corresponding entry in cellList).
%Possible actions include: 
%-Remove a field
%-Add a field and value
%-Rename a field
%-Reassign a value
%-Sync cel list
%-Delete the cell from cellList entirely

%CHANGE the project root if necessary (change for each machine)
physiologyRootPath = "C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/physiology";
addpath(genpath(physiologyRootPath)); %add all subfolders to matlab path

%CHANGE the file name of the list of all cells
cellListFName = "C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/batch/cellList_OKR.mat";

%% Select the cell that you would like to change
cd(physiologyRootPath)
[fname, fpath] = uigetfile('*.mat', 'Select the cell analysis structure you would like to edit');

load(fullfile(fpath, fname)) %load CellParameters structure for this cell to the workspace
load(cellListFName) %load cellList structure to the workspace
cellName = CellParameters.cellID;

%% Query operation that user would like to perform
operations = {'Remove a Field', 'Add a Field and Value',...
    'Rename a field', 'Reassign a Value', 'Sync cellList', 'Delete Cell Entirely'};

[selectionIndex] = listdlg('ListString', operations ,...
    'Name', 'Select an operation', 'ListSize', [300 300], 'SelectionMode', 'single');

% find where this cell is in cellListIndex
for i = 1 : size(cellList, 1)
    if strcmp(cellList{i, 1}, cellName)
        cellListIndex = i;
        break
    end
end


if selectionIndex < numel(operations) - 1 %last two operations don't work with this section
    
    allFields = fieldnamesr(CellParameters, 'full'); %generates cellArray of all fieldnames in structure
    allFields{end+1, 1} = 'CellParameters';
    %if removing a field, select the field to remove. If Adding a field,
    %select the field under which you'd like to add. To rename a field,
    %select the field to rename. To Reassign a value, select the field
    %underwhich the value is located. 
    fieldIndex = listdlg('ListString', allFields,...
        'Name', 'Select the field', 'ListSize', [600 600], 'SelectionMode', 'single');
    fieldName = allFields{fieldIndex};
    
    if strcmp(fieldName, 'CellParameters')
        %make fieldName a blank string if they want to edit the first level
        %of the structure. This makes indexing to it easier (so I don't
        %have to repeat the words "CellParameters" multiple times in the
        %code).
        fieldName = '';
    end
    
    
    if selectionIndex == 1
        %% Removing a field
        if numel(strfind(fieldName, '.')) > 0 %not a top level field
            
            %find index of last '.'
            reversedName = reverse(fieldName);
            lastDotIndex = numel(fieldName) - strfind(reversedName, '.') + 1;
            
            %grab names
            beforeDot = fieldName(1:lastDotIndex-1);
            afterDot = fieldName(lastDotIndex+1:end);
            
            %remove the field
            eval(strcat('CellParameters.', beforeDot,...
                '= rmfield(CellParameters.', beforeDot,',"', afterDot,'")'));
        
        else %top level field
           CellParameters = rmfield(CellParameters, fieldName);
        end
        
    elseif selectionIndex == 2 
        %% Add a field and its value below the selected field

        fieldAndValue = inputdlg({'Enter Field Name', 'Enter Value'}, 'Input');
        fieldToAdd = fieldAndValue{1};
        valueToAdd = fieldAndValue{2};
        
        if numel(fieldName) > 0 
            %Add a dot to the field name for proper indexing UNLESS the
            %user selected to edit in the top level of the structure (in
            %which case fieldName == '').
            fieldNameAdjusted = [fieldName, '.'];
        else
            fieldNameAdjusted = fieldName;
        end
                
        %if the value to add is not a string, change it to a double
        dType = questdlg('Select data type of value', 'Select Data Type', 'Num', 'String', 'String');
        switch dType
            case 'Num'
            accessString = strcat('CellParameters.', fieldNameAdjusted, fieldToAdd,...
            '=',valueToAdd);            
            
            otherwise   %value is a string (difference is "" around valueToAdd)
            accessString = strcat('CellParameters.', fieldNameAdjusted, fieldToAdd,...
            '= "',valueToAdd, '"');
        end
        
        eval(accessString);
        
    elseif selectionIndex == 3
        %% Rename a field      
        reversedName = reverse(fieldName);
        lastDotIndex = numel(fieldName) - strfind(reversedName, '.') + 1;

        %grab names
        beforeDot = fieldName(1:lastDotIndex-1);
        afterDot = fieldName(lastDotIndex+1:end);

        newFieldName = inputdlg('Enter new field name (must be valid matlab fieldname)', 'New Name', [1 50], {afterDot});
        newFieldName = newFieldName{1}; %grab string value       
        
        if numel(strfind(fieldName, '.')) > 0 %not a top level field

            %add a new field with that's a copy of the field it will
            %replace
            eval(strcat('CellParameters.', beforeDot, '.', newFieldName,...
                '= CellParameters.', fieldName));
            
            %Now remove the field with the old name
            eval(strcat('CellParameters.', beforeDot,...
                '= rmfield(CellParameters.', beforeDot,',"', afterDot,'")'));
            
        else %top levelfield
           CellParameters.newFieldName = CellParameters.fieldName;
           CellParameters = rmfield(CellParameters, fieldName);
        end
        
        
    elseif selectionIndex == 4
        %% Reassign a value
        newVal = inputdlg('Enter the new value to be assigned to your field');
        newVal = newVal{1};
        
        dType = questdlg('Select data type of value', 'Select Data Type', 'Num', 'String', 'String');
        
        switch dType
            case 'Num'
                accessString = strcat('CellParameters.', fieldName, '=', newVal);
            otherwise
                accessString = strcat('CellParameters.', fieldName, '= "', newVal, '"');
        end
        
        eval(accessString);
        
    end
    
end


%% executes for everything but deleting a cell entirely
if selectionIndex < numel(operations)
    cellList{cellListIndex, 2} = CellParameters; %syncs cellList with CellParameters
    cellList{cellListIndex, 2}.updateHTML = 1;
    
    %save CellParameters and cellList
    save(fullfile(fpath, fname), 'CellParameters');
    save(cellListFName, 'cellList');
    clear CellParameters cellList
end

%% executes if user wants to delete the cell entirely
if selectionIndex == numel(operations)
    cellList = [cellList(1:cellListIndex-1, :); cellList(cellListIndex+1:end, :)];
    save(cellListFName, 'cellList')
    warning(['PERMINENTLY deleted ', cellName, 'from cellList'])
end
    
    
    
    
    