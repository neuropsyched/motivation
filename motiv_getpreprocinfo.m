function [info] = motiv_getpreprocinfo(Nest)
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
%       dbs_subdir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    Nest = 'Trellis'; %Default Nest folder
end
info.projroot = '/Users/markrichardson/Documents/Projects/Motivation';
info.dataroot = '/Volumes/Nexus/Electrophysiology_Data/DBS_Intraop_Recordings/2014';
cd(info.projroot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp('Choose Patient Directory')
[info.patientdir] = uigetdir(info.dataroot,'Choose Patient Directory'); 
[~, filenames]    = dbs_subdir(fullfile(info.patientdir,Nest));
disp('All Trellis Filenames:')
disp(filenames')

info.reply = input('Which Blocks would you like to keep?','s'); % e.g. 234

tmp = cellfun(@(x) (x(1:strfind(filenames{1},'.')-1)), filenames,'UniformOutput',false); %remove extension
info.TrellisFilenames = horzcat({filenames{~cellfun(@isempty,strfind(tmp,info.reply(1)))}}); %first 3 cells
for i = 2:length(info.reply) %next cells by 3s
    info.TrellisFilenames = horzcat({info.TrellisFilenames{:}},{filenames{~cellfun(@isempty,strfind(tmp,info.reply(i)))}});
end
clc
end