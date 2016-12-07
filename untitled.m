%% import EMU data

ID=[1:64,129:160,257:288];
ch=ID(1:length(labels));
load('/Dropbox/Filters/filters.mat')
fs=30e3;

parfor i=1:length(ID)

    tic
    [~, raw, ~] = GetAnalogData([Path FileName '.ns5'], fs, ID(i));

    raw=bsxfun(@minus,raw,mean(raw,1));

    raw=unpowerline2(raw,fs,2);
    filt=filtfilt(low30k,raw);

    x=downsample(filt,3);

    filt=filtfilt(lp10k,x);
    x=downsample(filt,4);
    filtf(:,i)=filtfilt(hp2500,x);
    toc

end

filt=filtf;

eegplot(filt','srate',fs/12)
clearvars filtf
nfs=fs/12;
save('tmp.mat','-v7.3')

%% BIPOLAR REFERENCING
% All depth macroelectrodes (depth) should be bipolared. Micro electrodes should not be

% bipolared. 

elec={'LF','LGr','LPGr','LAS','LPS','LAD','LHD','LID' };
signal=zeros(size(filt));
temp2=[];
micro=[];
for e=1:length(elec)
    temp=find(~cellfun(@isempty,strfind(labels,elec{e})));
    if length(temp)==16
         micro=cat(1,micro,temp(9:end));
        temp=temp(1:8);
    end
    signal(:,temp(1:end-1))=bsxfun(@minus,filt(:,temp(1:end-1)),filt(:,temp(2:end)));
    temp2=cat(1,temp2,temp);
end
signal(:,micro)=filt(:,micro);

% Behnke Fried 
temp=find(~cellfun(@isempty,strfind(labels,'LBF')));
signal(:,temp)=bsxfun(@minus,filt(:,temp),filt(:,temp(length(temp))));
temp2=cat(1,temp2,temp);

%% COMMON AVERAGE REFERENCING

temp=setdiff(1:size(filt,2),temp2);
signal(:,temp)=bsxfun(@minus,filt(:,temp),mean(filt(:,temp),2));


eegplot(signal','srate',fs/12)
%% define the bad channels

badch=[];
filt(:,badch)=0;

%% ICA

[coeff,score,latent,tsquared,explained,~] = pca(filt(420*nfs:end,5:end));
plot(explained)
[weights,spheres]=runica(filt(420*nfs:end,5:end)','extended',1,'pca',5,'maxsteps',1e3);
eegplot(weights*spheres*filt(420*nfs:end,5:end)','srate',1200,'winlength',30,'dispchans',20)
[clean_data]=ahmad_icaproj(weights,spheres,filt(420*nfs:end,5:end)',[1 ]);
eegplot(ahmad_icaproj(weights,spheres,filt',[1 3]),'srate',1200,'winlength',10,'dispchans',20)

%% Psychometrics
%Gets Force
fs=1000;
[time1kHz, Force, IDs] = GetAnalogData([Path FileName '.ns2'], fs, [10241 10242]);

figure('windowstyle','docked'); plot(time1kHz, Force)

%% Two Trials
tstart= 1;
tend= round(580*nfs);
signal=signal(tstart:tend,:);
Force=Force(tstart:580*fs,:);
time1kHz=1/fs:1/fs:length(Force)/fs;

%%
tstart2= round(650*nfs);
signal3 = signal3(tstart2:end,:);
[time1kHz2, Force2, IDs] = GetAnalogData([Path FileName '.ns2'], fs, [10241 10242]);
Force2 = Force2(650*fs:end,:);


%%
[LeftResponseTimes,movementLt]  = MotivationTaskResponseTrig2(Force(:,1), 100, fs, 250, 100, 0.07);
[RightResponseTimes,movementRt]  = MotivationTaskResponseTrig2(Force(:,2), 100, fs, 250, 100, 0.07);


figure('windowstyle','docked'); plot(time1kHz, Force)
hold on; plot((1/fs)*LeftResponseTimes*[1 1], ylim, 'b');
hold on; plot((1/fs)*RightResponseTimes*[1 1], ylim, 'g');
hold on; plot((1/fs)*movementLt.deflection'*[1 1], ylim, 'r');
hold on; plot((1/fs)*movementRt.deflection'*[1 1], ylim, 'y');
scatter(movementLt.locs/fs,movementLt.peaks,'b','filled')
scatter(movementRt.locs/fs,movementRt.peaks,'r','d','filled')


%Timestamps for task parameters -- Just gotta know the order they appear
[EventTimes] = GetEventData([Path FileName '.nev']);
EventTimes = EventTimes(1:391,:);

[EventTimes2] = GetEventData([Path FileName '.nev']);
EventTimes2 = EventTimes2(392:end,:);
EventTimes2 = bsxfun(@minus,EventTimes2, 650);



figure('windowstyle','docked'); plot(time1kHz, Force)
hold on; plot(EventTimes*[1 1], ylim, 'k')
save('tmp.mat','-v7.3')

skips=7;

try
    celldata=BilateralDelayMotivationTaskIntraop;
    task=1;
catch
    celldata=BilateralDelayGripTaskIntraop;
    task=2;
end
%%
trialskips=1;
celldata(1:trialskips,:)=[];

%%
%Indexing event times
%Get BilateralDelayMotivationTaskIntraopData using import from txt

[CueStimulusTimes, CueStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 1, celldata);
[CommandStimulusTimes,~, CommandTrial, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 2, celldata);
[FeedbackStimulusTimes, ~, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 3, celldata);
[ITIStimulusTimes, ITIStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 4, celldata);
switch task  
    case 1
        RGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.Go));
        LGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.Go));
        RNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.NoGo));
        LNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.NoGo));
    case 2128
        RGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.Go));
        RGC=cat(2,RGC, intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.GoFast)));
        LGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.Go));
        LGC=cat(2,LGC2, intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.GoFast)));
        RNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Right,CommandTrial.NoGo));
        LNGC=intersect(CommandTrial.Correct,intersect(CommandTrial.Left,CommandTrial.NoGo));
end
RNGC=setdiff(RNGC,1);
LNGC=setdiff(LNGC,1);
RGC=setdiff(RGC,1);
LGC=setdiff(LGC,1);

for i=1:length(RGC)
    try
        RGCR(i,1)=RightResponseTimes(intersect(find(RightResponseTimes> fs*CommandStimulusTimes(RGC(i))),find( RightResponseTimes < fs*CommandStimulusTimes(RGC(i))+ 2*fs)));
    catch
        RGCR(i)=0;
       
    end
end
RGC=RGC(RGCR>0);
RGCR=RGCR(RGCR>0)/fs;
rind=cell2mat(arrayfun(@(x) intersect(find(movementRt.locs> fs*x),find( movementRt.locs < fs*x+ 2*fs)),RGCR,'UniformOutput',false));
movementRt.locs=movementRt.locs(rind);
movementRt.peaks=movementRt.peaks(rind);
movementRt.deflection=movementRt.deflection(rind);

for i=1:length(LGC)
    try
        LGCR(i,1)=     LeftResponseTimes(intersect(find(LeftResponseTimes> fs*CommandStimulusTimes(LGC(i))),find( LeftResponseTimes < fs*CommandStimulusTimes(LGC(i))+ 2*fs)));
    catch
        LGCR(i)=0;
  
    end
end
LGC=LGC(LGCR>0);
LGCR=LGCR(LGCR>0)/fs;
lind=cell2mat(arrayfun(@(x) intersect(find(movementLt.locs> fs*x),find( movementLt.locs < fs*x+ 2*fs)),LGCR,'UniformOutput',false));

movementLt.locs=movementLt.locs(lind);
movementLt.peaks=movementLt.peaks(lind);
movementLt.deflection=movementLt.deflection(lind);







[CueStimulusTimes2, CueStimulusDuration2, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 1, celldata2);
[CommandStimulusTimes2,~, CommandTrial2, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 2, celldata2);
[FeedbackStimulusTimes2, ~, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 3, celldata2);
[ITIStimulusTimes2, ITIStimulusDuration2, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 4, celldata2);
switch task
    case 1
        RGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.Go));
        LGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.Go));
        RNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.NoGo));
        LNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.NoGo));
    case 2
        RGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.Go));
        RGC2=cat(2,RGC2, intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.GoFast)));
        LGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.Go));
        LGC2=cat(2,LGC2, intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.GoFast)));
        RNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.NoGo));
        LNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.NoGo));
end
RNGC2=setdiff(RNGC2,1);
LNGC2=setdiff(LNGC2,1);
RGC2=setdiff(RGC2,1);
LGC2=setdiff(LGC2,1);

for i=1:length(RGC2)
    try
        RGCR2(i,1)=RightResponseTimes2(intersect(find(RightResponseTimes2> fs*CommandStimulusTimes2(RGC2(i))),find( RightResponseTimes2 < fs*CommandStimulusTimes2(RGC2(i))+ 2*fs)));
    catch
        RGCR2(i)=0;
       
    end
end
RGC2=RGC2(RGCR2>0);
RGCR2=RGCR2(RGCR2>0)/fs;
rind2=cell2mat(arrayfun(@(x) intersect(find(movementRt2.locs> fs*x),find( movementRt2.locs < fs*x+ 2*fs)),RGCR2,'UniformOutput',false));
movementRt2.locs=movementRt2.locs(rind2);
movementRt2.peaks=movementRt2.peaks(rind2);
movementRt2.deflection=movementRt2.deflection(rind2);

for i=1:length(LGC2)
    try
        LGCR2(i,1)=     LeftResponseTimes2(intersect(find(LeftResponseTimes2> fs*CommandStimulusTimes2(LGC2(i))),find( LeftResponseTimes2 < fs*CommandStimulusTimes2(LGC2(i))+ 2*fs)));
    catch
        LGCR2(i)=0;
  
    end
end
LGC2=LGC2(LGCR2>0);
LGCR2=LGCR2(LGCR2>0)/fs;
lind2=cell2mat(arrayfun(@(x) intersect(find(movementLt2.locs> fs*x),find( movementLt2.locs < fs*x+ 2*fs)),LGCR2,'UniformOutput',false));

movementLt2.locs=movementLt2.locs(lind2);
movementLt2.peaks=movementLt2.peaks(lind2);
movementLt2.deflection=movementLt2.deflection(lind2);
%% save
mkdir('/home/richardsonlab/Dropbox/STN_cortical connectivity/Epilepsy/AW093014')
save('/home/richardsonlab/Dropbox/STN_cortical connectivity/JM071415/JM08_RT.mat','-v7.3')

%%  Trial processing
for e=1:length(LGC)
%     TrialtypeL(e)=sum(ismember(CommandTrial.Win,LGC(e)))+2*sum(ismember(CommandTrial.Lose,LGC(e)))+3*sum(ismember(CommandTrial.Squeeze,LGC(e)));
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
[trialsL,badL]=check_trials(signal,ITIStimulusTimes(LGC-1),ITIStimulusTimes(LGC),ELeft,nfs);

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
    [pk loc]=max(Force(round(RGCR(e)*fs):round(RGCR(e)*fs)+fs,2));
    acc=(pk-Force(round(RGCR(e)*fs),2))/((loc-1)/fs);
    rxt=RGCR(128e)-CommandStimulusTimes(RGC(e));
    movement(e,:)=[rxt pk acc];
    ERight{e}=EEG;
    
end
[trialsR,badR]=check_trials(signal,ITIStimulusTimes(RGC-1),ITIStimulusTimes(RGC),ERight,nfs);

