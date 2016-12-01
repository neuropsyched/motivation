function [FileName] = motiv_getfilename(Path,ext)
% get filename from extension

files = dir(Path);
files = cat(2,{files(:).name});
fname = files(find(~cellfun(@isempty,strfind(files,ext))));

for i=1:size(fname,2)
    if ~isequal(fname{i}(1),'.')
        fname = fname{i};
    else
    end
end

idx = strfind(fname,ext)-1;
FileName = fname(1:idx);

end


% Notes:
% patientinfo = load(files(~isempyf

% files{find(~cellfun(@isempty,strfind(files,'data')))}

% if size(fnames,2)>1 && isequal(fnames{1},'.')
%     fnames = fnames{2}
% else