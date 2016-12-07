function [txtfiles,numtrials] = motiv_numtrials(filedir)
%% Import data from text file.
%% Initialize variables.
% clear
% clc
% PtId = 'RH021016';
% filedirectory = ['/home/richardsonlab/ari/Projects/' PtId '/Matlab/'];

if nargin<1;
    filedir = pwd;
end

[~,files] = dirdir(filedir);
txtfiles = files(~cellfun(@isempty,strfind(files,'.txt')));

txtfiles = txtfiles(listdlg('PromptString','Choose Files','SelectionMode','multiple','ListSize',[500,300],'ListString',txtfiles));

%
for n=1:length(txtfiles);

    % Open the text file.
    filename = txtfiles{n};
    fileID = fopen(filename,'r');      
    % Set textscan Options
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    delimiter = '\t';
    startRow = 3;
    % Read columns of data according to format string.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    % Close the text file.
    fclose(fileID);
    
    % Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    
    for col=[1,2,3,4,7,8,9,10,11,12]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(thousandsRegExp, ',', 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    
    % Split data into numeric and cell columns.
    rawNumericColumns = raw(:, [1,2,3,4,7,8,9,10,11,12]);
    rawCellColumns = raw(:, [5,6]);
    
    
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
    
    % Create output variable
    numtrials{n} = length(raw);
end
%% Clear temporary variables
% clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;