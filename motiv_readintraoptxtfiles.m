function BilateralDelayMotivationTaskIntraop = motiv_readintraoptxtfiles(filename)
%% Initialize variables.
% info = datastruct.info;
% try
%     [subdir,~] = dirdir(info.patientdir);
%     subdir = subdir(~cellfun(@isempty,strfind(lower(subdir),'matlab')));
%     cd(char(fullfile(info.patientdir,subdir)))
% catch
% end
% info.MTdir = char(fullfile(info.patientdir,subdir));
% [~,files] = dirdir(info.MTdir);
% files = files(~cellfun(@isempty,strfind(files,'intraop')));
% files = files(~cellfun(@isempty,strfind(files,'.txt')));
%
% f=1;
% files = files(listdlg('PromptString',['Choose Text Files: for ' datastruct.FilesNames{f}],...
%     'SelectionMode','multiple','ListSize',[450 300],'ListString',files));
% info.MT_files{f} = files;
%
% for i = 1:length(files)
if ~ischar(filename)
    filename = char(filename);
end
delimiter = '\t';
startRow = 3;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
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

rawNumericColumns = raw(:, [1,2,3,4,7,8,9,10,11,12]);
rawCellColumns = raw(:, [5,6]);

BilateralDelayMotivationTaskIntraop = raw;
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;
end