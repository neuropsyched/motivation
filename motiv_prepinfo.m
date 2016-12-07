function [info] = motiv_prepinfo(PtId)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Get all Trellis File Names folder
%
% Inputs:
%           Nest - folders ./patientdirectory containing .NSX files (ONLY)
% Outputs:
%           info.projroot: USER entered path to project root
%           info.dataroot: USER entered path to root for all data
%           info.patinetstr: USER enetered path
%
% Dependent Functions:
%       dirdir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.projroot = fileparts(which('motiv_prepinfo.m'));
info.dataroot = '/home/richardsonlab/ari/Projects';
% cd(info.projroot);

if nargin<1
    disp('Choose Patient Directory')
    [subs,~] = dirdir(info.dataroot);
    info.patientname = subs{listdlg('PromptString','Choose Patient','SelectionMode','single','ListString',subs)};
elseif nargin==1
    if ~ischar(PtId)
        error('Patient ID must be a string')
    end
    info.patientname = PtId;
end
if isempty(info.patientname)
    error('Please select a patient from the list or enter one as an argument')
end

info.patientdir = fullfile(info.dataroot,info.patientname);

% Check if directory exists
if ~exist(info.patientdir,'dir')
    error('Nonexistent Patient Directory')
end
% cd(info.patientdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Pass Enter All Information and Save motiv_info.mat
if ~exist(fullfile(info.patientdir,'motiv_info.mat'),'file')
clc
[foldernames,filenames]    = dirdir(info.patientdir);
info.trellispath = (info.patientdir);
if sum(~cellfun(@isempty,strfind(filenames,{'.ns2'})))==0
    idx = find(~cellfun(@isempty,strfind(lower(foldernames),'trellis'))==1);
    [~,filenames]    = dirdir(fullfile(info.patientdir,char(foldernames(idx))));
    info.trellispath = fullfile(info.patientdir,char(foldernames(idx)));
end
clear idx
disp('All Trellis Filenames:')
disp(filenames')

% Filesets to keep
info.reply = strsplit(input('Enter the (Comma Separated) Blocks that you would like to keep? ','s'),','); % e.g. 2,3,4
% info.reply = strsplit(cell2mat(inputdlg('Enter the blocks that you would like to keep')),','); % :)


% get temporary filenames without prepend or extension
for i = 1:length(filenames)
    idx = strfind(filenames,'00');
    tmp{i} = filenames{i}(idx{i}(1):end-4);
end

% get TrellisFilenames by group
info.TrellisFilenames = {};
for i = 1:length(info.reply) %next cells by 3s
    info.TrellisFilenames = horzcat({info.TrellisFilenames{:}},{filenames{~cellfun(@isempty,strfind(tmp,info.reply{i}))}});
end

clc
save(fullfile(info.patientdir,'motiv_info.mat'),'info')
elseif exist(fullfile(info.patientdir,'motiv_info.mat'),'file')
    load(fullfile(info.patientdir,'motiv_info.mat'))
end

end

% if ~cellfun(@isempty,strfind(foldernames,'Trellis'))
% ~cellfun(@isempty,strfind(filenames,{'.nev'})) + ...
%     ~cellfun(@isempty,strfind(filenames,{'.ns2'})) + ...
%         ~cellfun(@isempty,strfind(filenames,{'.ns5'})) 

% CELLFUN
% tmp = cellfun(@(x) (x(1:strfind(filenames{1},'.')-1)), filenames,'UniformOutput',false); %remove extension
% info.TrellisFilenames = horzcat({filenames{~cellfun(@isempty,strfind(tmp,info.reply{1}))}}); %first 3 cells
