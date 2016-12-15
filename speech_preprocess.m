%%                   Initialize Variables
%%%%%%%%%%%%%%%%
%% import LFP data
clear; clc
info = motiv_prepinfo;
[xlpath,xlname,xlext] = fileparts('~/Documents/Projects/ECoG Electrode Keys.xlsx');
xlsheet = info.patientname(3:end);
try
    [~, ~, Labels] = xlsread(fullfile(xlpath,[xlname,xlext]),xlsheet);
    info.labels = Labels(4:129,2);
    info.labels(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),info.labels)) = {''};
catch
end
clear xl*
% IDs=[1:64,129:160,257:288];
% info.ch=IDs(:);
info.fs=1000;
info.nev = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.nev')));
info.ns2 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns2')));
info.ns5 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns5')));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[datastruct] = motiv_getfiledata(info); % block=3; ID=[]; 

%% define the bad channels
ch=datastruct.ChannelIDs;
f=1; % file index
filt=datastruct.filtdata{f};
V=var(filt);
figure;
scatter(1:length(ch),V)
badch=[29 30 31 32];
%% referencing
%CAR
fs=1200;
% signal=bsxfun(@minus,filt,mean(filt(:,setdiff(1:length(ch),badch)),2));
% eegplot(eegfilt(signal',fs,1,0),'srate',fs)
signal=bsxfun(@minus,filt(:,setdiff(1:length(ch),badch)),mean(filt(:,setdiff(1:length(ch),badch)),2));
eegplot(signal','srate',fs,'title',datastruct.FileNames{f})
datastruct.signal{1}=signal;

%% import NO data
mat_files = dir('*.mat');
input={'CRAW', 'CMacro_LFP'};
output={'Raw','LFP'};
[Rec]=ReadNeuroOmega_AA(input, output,mat_files);

fName = fullfile(datastruct.info.patientdir,[datastruct.info.patientname '_NO_Rec_raw.mat']);
save(fName,'Rec','-v7.3')
clear input output mat_files fName

%% preprocess NO data 
Filters = load(fullfile(datastruct.info.projroot,'Filters/44kfilters.mat'));
lp1200_44k = Filters.lp1200_44k;
lp500_4k = Filters.lp500_4k;
filt_micro = cell(size(Rec)); filt_macro = filt_micro;
filt_micronspk=[];

% Raw data
parpool(4)
parfor i =1: length(Rec)
    figure
    tmp=bsxfun(@minus,Rec(i).Raw.ts,mean(Rec(i).Raw.ts,2))';
    plot_fft(Rec(i).Raw.ts(1,:),Rec(i).Raw.sr);
    hold on
    tmp=unpowerline(tmp,Rec(i).Raw.sr,1);
    plot_fft(tmp(1,:),Rec(i).Raw.sr);
    
    tmp=filtfilt(lp1200_44k,tmp);
    
    ds=downsample(tmp,11);
    tmp=filtfilt(lp500_4k,ds);
    plot_fft(tmp(1,:),Rec(i).Raw.sr/11);
    
    ds=downsample(tmp,2);
    filt_micro{i}=ds;
end
srmic=2e3;

% Macro LFP data
for i =1: length(Rec)
    tmp=bsxfun(@minus,Rec(i).LFP.ts,mean(Rec(i).LFP.ts,2))';
    tmp=unpowerline(tmp,Rec(i).LFP.sr,1);
    filt_macro{i}=resample(tmp,1200,1375)';
end

%% Align ECOG and NO 
file = 1;%left

%
fileName = fullfile(datastruct.info.trellispath,char(datastruct.info.nev{file}));
[EventTimesTrellis] = GetEventData(fileName); % trellis time stamps 
clear fileName

%
figure; plot(EventTimesTrellis*[1 1],ylim)
nfs = 1200;
filt = datastruct.filtdata{file};
% visually identify T0 / T1 in Trellis as the start / end of segment of interest
T0 = [3200 2600 2275 1550];
for idx=1:length(Rec)
    ET{idx} = sort([Rec(idx).DigUp Rec(idx).DigDown]);  % NO event times
    
    Event0 = EventTimesTrellis(find(EventTimesTrellis>T0(idx),1,'first'));
    tstart(idx)= Event0 - ET{1,idx}(1);
    try
        tend(idx) = tstart(idx) + length(Rec(idx).Raw.ts)/Rec(idx).Raw.sr; % 5 secs after last time stamp
        E{idx} = filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);
    catch
        tend(idx) = tstart(idx) + ET{1,idx}(end)+10; % 5 secs after last time stamp
        E{idx} = filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);
    end
    
end

%%  save into task files
fileName=[datastruct.info.nev{file}(9:end-4) '_raw.mat'] ;
save(fileName,'E','filt_macro','fileName','nfs','srmic',...
    'ET','T0','filt','EventTimesTrellis','Rec','datastruct','-v7.3')

