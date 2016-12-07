function [filestr] = motiv_getpatientlist(root)

%small function to get patientlist
try
    filestr = load(fullfile(root,'patientlist.mat'));
catch
    fName = '/Users/markrichardson/Google Drive/Richardson Lab/Databases/MotivationTask_Database_103116.xlsx';
    [~, ~, filestr] = xlsread(fName,'DBSLocsMotiv');
    filestr = filestr(2:24,1);
    filestr(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),filestr)) = {''};
end