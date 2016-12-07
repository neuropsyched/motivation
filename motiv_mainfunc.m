% function motiv_mainfunc
%% motiv_open
clear
clc

% Find Root Directory
% Add to varinfo (variable information)
[globalinfo.root,~,~] = fileparts(which('motiv_mainfunc.m'));
cd(globalinfo.root)

try
    load(fullfile(root,'motivinfo.mat'))
catch
    load(which('motivinfo.mat'))
end

% other options
globalinfo.recordstr = {motivinfo.pID};

globalinfo.fileHbipolar = {'_Cmd','_ITI','_Rsp'}; % file handle for Motivation_Bipolar
globalinfo.fileHraw = {'_cavref','_noref','_raw','_trial'}; % file handle for Motivation_Raw



%% RUN fbs Functions
clear matObj*
matObj = matfile('Motivation_Bipolar','writable',true);
[objOut,psdOut,sList] = fbs_trialbypsd(matObj);
% clear matObj
% save('psdOutput.mat','psdOut')
% [m,n] = fbs_trialbypsd(matObj)

fbs_plotpsd(objOut,psdOut,sList)








%% Load Bipolar Variables of Interest & Set Local Variable Environment
% USER INPUT
mType = 'Bipolar'; % mType is either 'Raw' or 'Bipolar'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE
if ~strcmp(mType,'Raw') && ~strcmp(mType,'Bipolar') && isequal(mType,2)
    error('enter string for file type (mType): ''Raw'' or ''Bipolar''')
end
% if isequal(mType,2) % then run both in a future version
mName = fullfile(globalinfo.root,['Motivation_' mType '.mat']); %Name of Matfile to search
if ~exist(mName)
    error(['missing file - ' mName])
end
fileH = eval(['globalinfo.fileH' lower(mType)]); %local fileHandle
% Select Recordings
[sList,~] = listdlg('ListString',globalinfo.recordstr,...
    'PromptString','Select recordings to analyze','SelectionMode','multiple');
fName = cell(size(sList));
for i = 1:length(sList) %number of recordings selected
    for j=1:length(fileH);
        fName{i} = [globalinfo.recordstr{sList(i)} fileH{j}];
        eval([char(fName{i}) ' = load(''' mName ''',''' char(fName{i})  ''');'])
        while iscell(eval(char(fName{i}))) || isfield(eval(char(fName{i})),char(fName{i}))
            if iscell(eval(char(fName{i})))
                hold = eval(char(fName{i}));
                eval([char(fName{i}) ' = hold{1};'])
                clear hold
            elseif isfield(eval(char(fName{i})),char(fName{i}))
                hold = eval(char(fName{i}));
                eval([char(fName{i}) ' = {hold};'])
                clear hold
            end
        end
                    
    end
end
% save local variables
local.mType = mType;
local.mName = mName;
local.fileH = fileH;
local.sList = sList;
local.fName = fName;
clear ans i j mName fName fileHandle sList mType





%% PLOT PSD using Local Variables
% plot(squeeze(LM092314_1_ITI(:,12,:))) %[
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
if ~strcmp(local.mType,'Bipolar')
    error('Cannot use Raw data here')
end
ntrial = 12; % plot trial 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE
i=1;
% data = eval([globalinfo.recordstr{local.sList(i)} local.fileH{i}]); % see i
srate = motivinfo(local.sList).fs; %find srate from struct, assume its same in both recordings
time = (0:size(data,1)-1)/srate; %time is in the first dimension

% Remove NANs
if ~isempty(isnan(data));
    nanidx = find(isnan(data));
    data(nanidx)=0;
end

figure;
plot(time,squeeze(data(:,ntrial,:)))
set(gca,'xlim',[0 size(data,1)/srate])
    xlabel('Time(s)')
    ylabel('Frequency(Hz)')
    title('Time Domain')
    grid('on')

figure;
pwelch(squeeze(data(:,ntrial,:))) %
% clear hold_data
% Are there NANs in the data?: y = ~isempty(isnan(hold_data));
clear h i n

%% load cavref noref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
id = 'LM'; % identifiable string to load patient data
h = 1;  % if 2 datasets, choose the first 1
ntrial = 12; % plot trial 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE
n = find(~cellfun(@isempty,strfind(filestr,id))==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h = 1:h
    for i=1:length(fileH2);
        eval(['load(''' filestr{n(h)} fileH2{i} ''')'])
    end
end
clear h i n

data = eval([filestr{n(h)} fileH{i}]); % see i
[srate,~] = pinfo(indx.pinfo(n)).fs; %find srate from struct, assume its same in both recordings
time = (0:size(data,1)-1)/srate; %time is in the first dimension

%% Load Raw Variables of Interest
% matObj = matfile('Movement_Bipolar.mat','writable',true)
% filestr(find(~cellfun(@isempty,strfind(filestr,'LM'))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
id = 'LM'; % identifiable string to load patient data
h = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE
n = find(~cellfun(@isempty,strfind(filestr,id))==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bName = 'Motivation_Data';
for h = 1:h
    for i=1:length(fileH2);
        eval([filestr{n(h)} fileH2{i} ' = load(''' filestr{n(h)} fileH2{i} ''');'])
    end
end
clear h i n
% clear *1* *2*
