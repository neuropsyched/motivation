function [var1,var2,var3] = motiv_loadsets(globalinfo,mType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Small function to load data from object
% pass globalinfo from motiv_mainfunc.m
% USER INPUT
% id = 'LM'; % identifiable string to load patient data
% recording = 1;
% mName = 'Motivation_Data.mat']; %Name of Matfile to search
% mType is either 'Raw' or 'Bipolar'
% =========================================================================
% Ari Kappel, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(mType,'Raw') || ~strcmp(mType,'Bipolar')
    error('enter string for file type (mType): ''Raw'' or ''Bipolar''')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE
mName = fullfile(globalinfo.root,['Motivation_' mType '.mat']); %Name of Matfile to search
if ~exist(mName)
    error(['missing file - ' mName])
end
fileHandle = eval(['globalinfo.fileH' lower(mType)]);
%% Select Recordings
[sList,~] = listdlg('ListString',globalinfo.recordstr,...
    'PromptString','Select recordings to analyze','SelectionMode','multiple');
for i = 1:length(sList) %number of recordings selected
    for j=1:length(fileHandle);
        fName = [globalinfo.recordstr{sList(i)} fileHandle{j}];
        eval([fName ' = load(''' mName ''',''' fName  ''');'])
    end
end
end



% % search mfind index of patient identifier (id)
% index = find(~cellfun(@isempty,strfind(globalinfo.recordstr,id))==1); 
% %%%%%%%%%%%%%%%%%%%
% for h = 1:nrecords %number of recordings selected
%     for i=1:length(motivinfo.varinfo.fileH);
%       fName = [globalinfo.filestr{index(h)} globalinfo.fileH{i}];
%       eval(['load(''' mName ''',''' fName  ''');'])
%     end
% end
% end