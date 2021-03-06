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

save('datastruct.mat','datastruct','-v7.3')
%% define the bad channels
load('datastruct.mat','datastruct','-v7.3')
datastruct.nfs = 1200;
ch = datastruct.ChannelIDs;
for f = 1:length(datastruct.FileNames)
    % f=1; % file index
    filt=datastruct.filtdata{f};
    V=var(filt);
    figure;
    scatter(1:length(ch),V)
    keyboard
    datastruct.badch{f}=[29 30 31 32];
end
clear V f filt ch ans
%% referencing
for f = 1:length(datastruct.FileNames);
    ch = datastruct.ChannelIDs;
    nfs = datastruct.nfs;
    badch = datastruct.badch{f};
    filt=datastruct.filtdata{f};
    % signal=bsxfun(@minus,filt,mean(filt(:,setdiff(1:length(ch),badch)),2));
    % eegplot(eegfilt(signal',fs,1,0),'srate',fs)
    signal=bsxfun(@minus,filt(:,setdiff(1:length(ch),badch)),mean(filt(:,setdiff(1:length(ch),badch)),2));
    eegplot(signal','srate',nfs,'title',datastruct.FileNames{f});
    datastruct.signal{f}=signal;
    keyboard
end
%% import NO data
input={'CRAW', 'CMacro_LFP'};
output={'Raw','LFP'};
for f = 1:length(datastruct.FileNames)
    sides = {'Left' 'Right'};
    cd([datastruct.info.patientdir '/NeurOmega/' sides{f}])
    mat_files = dir('*.mat');
    [Rec] = ReadNeuroOmega_AA(input, output, mat_files);
    Rec = Rec(end:-1:1);
    datastruct.Rec{f} = Rec;
    clear Rec mat_files
end

% fName = fullfile(datastruct.info.patientdir,[datastruct.info.patientname '_NO_Rec_raw.mat']);
% save('datastruct2.mat','datastruct','-v7.3')
clear input output mat_files fName

%% preprocess NO data 
Filters = load(fullfile(datastruct.info.projroot,'Filters/44kfilters.mat'));
lp1200_44k = Filters.lp1200_44k;
lp500_4k = Filters.lp500_4k;
nfs = datastruct.nfs;
% Raw data
parpool(8)
for f = 1:length(datastruct.FileNames);
    Rec = datastruct.Rec{f};
    filt_raw = cell(size(Rec)); filt_lfp = cell(size(Rec)); 
    srmic=2e3;
        
    % Micro Raw
    parfor i =1:length(Rec)
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
        filt_raw{i}=ds;
    end
    datastruct.filt_Raw{f}=filt_raw;
end
for f = 1:length(datastruct.FileNames);
    Rec = datastruct.Rec{f};
    filt_raw = cell(size(Rec)); filt_lfp = cell(size(Rec)); 
    srmic=2e3;
    % Macro LFP data
    for i =1: length(Rec)
        tmp=bsxfun(@minus,Rec(i).LFP.ts,mean(Rec(i).LFP.ts,2))';
        tmp=unpowerline(tmp,Rec(i).LFP.sr,1);
        filt_lfp{i}=resample(tmp,nfs,1375)';
    end
    datastruct.filt_LFP{f}=filt_lfp;
    datastruct.srmic{f}=srmic;

end
save('datastruct2.mat','datastruct','-v7.3')
%% Align ECOG and NO 
% load('datastruct2.mat')
nfs = datastruct.nfs;
fs = 1000;
ask = 1;
for f = 1:length(datastruct.Rec)
    Rec = datastruct.Rec{f}; %Previously flipped to decreasing depths to match file order
    filt = datastruct.filtdata{f};
    force = datastruct.forcedata{f};
    %
    fileName = fullfile(datastruct.info.trellispath,char(datastruct.info.nev{f}));
    [EventTimesTrellis] = GetEventData(fileName); % trellis time stamps     
    
    switch ask
        case 1
        h = figure; plot(EventTimesTrellis*[1 1],ylim)
        disp(['Displaying EventTimes for datfile' datastruct.FileNames{f}])
        title(['datfile' num2str(datastruct.FileNames{f})])
        waitfor(h) % pause while figure is open

        disp('Choose start and stop times')
        % visually identify T0 / T1 in Trellis as the start / end of segment of interest
        T0 = str2double(strsplit(cell2mat(inputdlg('Enter start times (space delimited)')),' '));
        
        case 0
            switch f
                case 1
                   T0 = [1700,2300,2700,3620];
                case 2
                    T0 = [280,740,1530];
            end
    end

    for idx=1:length(Rec)
        EventTimesNO{idx} = sort([Rec(idx).DigUp Rec(idx).DigDown])';  % NO event times
        % h=figure('windowstyle','docked'); plot(EventTimesNO{idx}*[1 1],ylim); title('EventTimesNO');
        % h=figure('windowstyle','docked'); plot(EventTimesTrellis*[1 1],ylim); title('EventTimesTrellis');
        
        Event0 = EventTimesTrellis(find(EventTimesTrellis>T0(idx),1,'first'));
        tstart(idx)= Event0 - EventTimesNO{idx}(1);
        try
            tend(idx) = tstart(idx) + length(Rec(idx).Raw.ts)/Rec(idx).Raw.sr; % 5 secs after last time stamp
            E{idx} = filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);
            F{idx} = force(round(tstart(idx)*fs):round(tend(idx)*fs),:);
        catch
            tend(idx) = tstart(idx) + EventTimesNO{idx}(end)+10; % 5 secs after last time stamp
            E{idx} = filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);
            F{idx} = force(round(tstart(idx)*fs):round(tend(idx)*fs),:);
        end
    end
    datastruct.EventTimesTrellis{f} = EventTimesTrellis;
    datastruct.EventTimesNO{f} = EventTimesNO;
    datastruct.E{f} = E;
    datastruct.F{f} = F;
    clear EventTimesTrellis EventTimesNO E F 
end
clear ask f filt force idx tend tstart T0 Event0 Rec
% save('datastruct3.mat','datastruct','-v7.3')
% filtered_data = datastruct.E;
% save('filtered_data_aligned.mat','filtered_data')

%% Align responses
% load('/home/ari/Documents/Projects/Motivation/RH021016/datastruct3.mat');
info = datastruct.info;
fs=1000;
nfs = datastruct.nfs;
[subdir,~] = dirdir(info.patientdir);
subdir = subdir(~cellfun(@isempty,strfind(lower(subdir),'matlab')));
info.MTdir = char(fullfile(info.patientdir,subdir));

if ~exist(info.MTdir,'dir')
    error(['Nonexistent Directory: Tried to access: ' info.MTdir])
end
[~,files] = dirdir(info.MTdir);
files = files(~cellfun(@isempty,strfind(files,'intraop')));
files = files(~cellfun(@isempty,strfind(files,'.txt')));
for f=1:length(datastruct.FileNames);
    MT_files = files(listdlg('PromptString',['Choose Text Files for: ' datastruct.FileNames{f}],...
        'SelectionMode','multiple','ListSize',[450 300],'ListString',files));
    info.MT_files{f} = MT_files;
end

f=1; %length(datastruct.FileNames) %file
t=1; %length(info.MT_files{f}) %num_task
celldata = motiv_readintraoptxtfiles(fullfile(info.MTdir,info.MT_files{f}(t))); 

% isempty(find(cellfun(@isempty,find(celldata(:,7),'ABORT'))==1))
% try
%     celldata=BilateralDelayMotivationTaskIntraop;
% end
trialskips = 1; %no baseline for the first trial
celldata(trialskips,:)=[];
%Indexing event times
%Get BilateralDelayMotivationTaskIntraopData using import from txt
skips = 7;
EventTimes=datastruct.EventTimesNO{f}{t}; % figure; plot(EventTimes*[1 1],ylim)
[CueStimulusTimes, CueStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 1, celldata);
[CommandStimulusTimes,~, CommandTrial, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 2, celldata);
[FeedbackStimulusTimes, ~, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 3, celldata);
[ITIStimulusTimes, ITIStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 4, celldata);
task=1;
switch task
    case 1
        RGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.Go));
        LGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.Go));
        RNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.NoGo));
        LNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.NoGo));
    case 2
        RGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.Go));
        RGC=cat(2,RGC, intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.GoFast)));
        LGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.Go));
        LGC=cat(2,LGC, intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.GoFast)));
        RNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.NoGo));
        LNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.NoGo));
end
RNGC=setdiff(RNGC,1);
LNGC=setdiff(LNGC,1);
RGC=setdiff(RGC,1);
LGC=setdiff(LGC,1);

for i=1:length(RGC)
    try
        RGCR(i,1)=     RightResponseTimes(intersect(find(RightResponseTimes> fs*CommandStimulusTimes(RGC(i))),find( RightResponseTimes < fs*CommandStimulusTimes(RGC(i))+ 2*fs)));
    catch
        RGCR(i)=0;
    end
end
RGC=RGC(RGCR>0);
RGCR=RGCR(RGCR>0)./fs;
for i=1:length(LGC)
    try
        LGCR(i,1)=     LeftResponseTimes(intersect(find(LeftResponseTimes> fs*CommandStimulusTimes(LGC(i))),find( LeftResponseTimes < fs*CommandStimulusTimes(LGC(i))+ 2*fs)));
    catch
        LGCR(i)=0;
    end
end
LGC=LGC(LGCR>0);
LGCR=LGCR(LGCR>0)./fs;

%% spectral measures
amp=abs(hilbert(eegfilt(filt',fs,70,200)'));
ampC=abs(hilbert(eegfilt(CAR',fs,70,200)'));
ampB=abs(hilbert(eegfilt(refbi',fs,80,150)'));

%%  Trial processing
f=1;
t=1;
filtered = datastruct.E{f}{t};
for e=1:length(LGC)
    % TrialtypeL(e)=sum(ismember(CommandTrial.Win,LGC(e)))+2*sum(ismember(CommandTrial.Lose,LGC(e)))+3*sum(ismember(CommandTrial.Squeeze,LGC(e)));
    EEG.event(1).type='Cue';
    EEG.event(1).latency=fs*(CueStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    EEG.event(2).type='Command';
    EEG.event(2).latency=fs*(CommandStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    EEG.event(3).type='Squeeze';
    EEG.event(3).latency=fs*(LGCR(e)-ITIStimulusTimes(LGC(e)-1));
    EEG.event(4).type='Feedback';
    EEG.event(4).latency=fs*(FeedbackStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    ELeft{e}=EEG;
end
[trialsL,CTIL]=check_trials(filtered,ITIStimulusTimes(LGC-1),ITIStimulusTimes(LGC),ELeft,fs);

for e=1:length(RGC)
    TrialtypeR(e)=sum(ismember(CommandTrial.Win,RGC(e)))+2*sum(ismember(CommandTrial.Lose,RGC(e)))+3*sum(ismember(CommandTrial.Squeeze,RGC(e)));
    EEG.event(1).type='Cue';
    EEG.event(1).latency=fs*(CueStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    EEG.event(2).type='Command';
    EEG.event(2).latency=fs*(CommandStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    EEG.event(3).type='Squeeze';
    EEG.event(3).latency=fs*(RGCR(e)-ITIStimulusTimes(RGC(e)-1));
    EEG.event(4).type='Feedback';
    EEG.event(4).latency=fs*(FeedbackStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    [pk loc]=max(Force(round(RGCR(e)*fs):round(RGCR(e)*fs+fs),2));
    acc=(pk-Force(round(RGCR(e)*fs),2))/((loc-1)/fs);
    rxt=RGCR(e)-CommandStimulusTimes(RGC(e));
    movement(e,:)=[rxt pk acc];
    ERight{e}=EEG;
    
end
[trialsR,badRT]=check_trials(reshape(smooth(ampC,50),[],length(labels)),ITIStimulusTimes(RGC-1),ITIStimulusTimes(RGC),ERight,fs);

rxn=mean(RGCR-CommandStimulusTimes(RGC)');
for e=1:length(RNGC)
    TrialtypeR(e)=sum(ismember(CommandTrial.Win,RNGC(e)))+2*sum(ismember(CommandTrial.Lose,RNGC(e)))+3*sum(ismember(CommandTrial.Squeeze,RNGC(e)));
    EEG.event(1).type='Cue';
    EEG.event(1).latency=fs*(CueStimulusTimes(RNGC(e))-ITIStimulusTimes(RNGC(e)-1));
    EEG.event(2).type='Command';
    EEG.event(2).latency=fs*(CommandStimulusTimes(RNGC(e))-ITIStimulusTimes(RNGC(e)-1));
    EEG.event(3).type='Squeeze';
    EEG.event(3).latency=fs*(EEG.event(2).latency/fs + rxn);
    EEG.event(4).type='Feedback';
    EEG.event(4).latency=fs*(FeedbackStimulusTimes(RNGC(e))-ITIStimulusTimes(RNGC(e)-1));
    ERightNG{e}=EEG;
end
[trialsRNG,CTIRNG]=check_trials(filtered,ITIStimulusTimes(RNGC-1),ITIStimulusTimes(RNGC),ERightNG,fs);
