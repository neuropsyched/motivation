%% import LFP data
[xlname,xlpath,~] = uigetfile({'*.xlsx'},'Select Labels Spreadsheet');
A=importdata(fullfile(xlpath,xlname));
[labels,ch]=xl_labels(A);
IDs=[1:64,129:160,257:288];
ch=IDs(97:97+63);
fs=1000;
[~, raw, IDs] = GetAnalogData([Path FileName '.ns2'], fs, ch);

%% data preprocessing
raw=bsxfun(@minus,raw,mean(raw,1));
% notched=unpowerline2(raw,fs,0.5);
filt=filtfilt(low,raw);
filt=filtfilt(high,filt);
eegplot(filt','srate',fs)
clearvars x raw notched
%% define the bad channels
badch=[51:58 19:26];

filt(:,badch)=[];
labels(badch)=[];
%% referencing
% bipolar
tf = isstrprop(labels, 'alpha');
for i=1:length(labels)
    elec{i}=labels{i,1}(tf{i,1});
end
elec=unique(elec);
refbi=zeros(size(filt));
for e=1:length(elec)
    temp=find(~cellfun(@isempty,strfind(labels,elec{e})));
    refbi(:,temp(1:end-1))=bsxfun(@minus,filt(:,temp(1:end-1)),filt(:,temp(2:end)));
end
eegplot(refbi','srate',fs)

for i=[1:32]
    refavg(:,i)=filt(:,i)-filtered(:,11);
end

for i=1:4
    refavg(:,i)=filtered(:,i)-filtered(:,10);
end

%CAR
CAR=bsxfun(@minus,filt,mean(filt,2));

CAR(:,1:94)=bsxfun(@minus,filt(:,1:94),mean(filt(:,1:94),2));
CAR(:,95:102)=bsxfun(@minus,filt(:,95:102),mean(filt(:,95:102),2));
CAR(:,103:108)=bsxfun(@minus,filt(:,103:108),mean(filt(:,103:108),2));
eegplot(CAR','srate',fs)
%% ICA
[~,~,~,~,explained,~] = pca(filt);
plot(explained)
[weights,spheres]=runica(filt','extended',1,'pca',20,'maxsteps',1e3);
[clean_data]=ahmad_icaproj(weights,spheres,filt',1);

%%
%Gets Force
fs=1000;
[time1kHz, Force, IDs] = GetAnalogData([Path FileName '.ns2'], fs, [10241 10242]);
figure; plot(time1kHz, Force)

getForceDataAA
hold on; plot((1/fs)*LeftResponseTimes*[1 1], ylim, 'b');
hold on; plot((1/fs)*RightResponseTimes*[1 1], ylim, 'g');
%Timestamps for task parameters -- Just gotta know the order they appear
[EventTimes] = GetEventData([Path FileName '.nev']);
hold on; plot(EventTimes*[1 1], ylim, 'k')
%% Align responses
trialskips=1;
skips=7;
try
    celldata=BilateralDelayMotivationTaskIntraop;
    task=1;
catch
    celldata=BilateralDelayGripTaskIntraop;
    task=2;
end
celldata(trialskips,:)=[];
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
%% Reject bad trials
trialsR(badRT)=[];
ERight(badRT)=[];
TrialtypeR(badRT)=[];
movement(badRT,:)=[];
%% find task reactive channels
Events=RGCR;
trial=RGC;
prestim=0.2;
poststim=1.5;
if sum(badRT)~=0
trial(badRT)=[];
Events(badRT)=[];
end
baseline=ITIStimulusTimes(trial-1);
baseline_dur=min(ITIStimulusDuration(trial));
time=-prestim:1/fs:poststim;
[package]=packagedata(reshape(smooth(ampC,50),[],length(labels)),fs,Events,prestim,poststim,baseline,baseline_dur);

[R,GR,W,L,S]= task_rxt_ch(package,TrialTypeR,100);
for i=1:size(GR,1)
    figure
    plot(time,mean(W(GR{i,1},:,:),3))
    hold on
    plot(time,mean(L(GR{i,1},:,:),3),'r')
    plot(time,mean(S(GR{i,1},:,:),3),'y')
    for t=1:size(GR{i,2},2)
        plot(time(GR{i,2}{t}),ones(length(GR{i,2}{t}),1))
    end
    title(labels(GR{i,1}))
    
end

[Corr]= gamma_corr(package,[],[0.2*fs (0.2*fs)+500],movement);
[Corrsp]= gamma_corr(package,TrialType,[0.2*fs (0.2*fs)+500],movement);


%% phase encoding
% define trials trialtype prestim poststim prebase Einfo EM
TrialType=TrialtypeR;
[complex]=packagedata(hilbert(eegfilt(CAR(:,elec)',fs,4,8)'),fs,Events,prestim,poststim,baseline,baseline_dur);
complex.trials=permute(complex.trials,[2 3 1]);
complex.baseline=permute(complex.baseline,[2 3 1]);
Condition=ones(length(TrialType),1);
Condition(find(TrialType==2))=2;
Condition(find(TrialType==3))=3;
[Epairs,pval,rho]=phase_encoding(angle(complex.trials),Condition);
[Epairsb,pvalb,rhob]=phase_encoding(angle(complex.baseline),Condition);
% normalize and group by region
normrho=bsxfun(@minus,rho,mean(rhob,1));
intraPFC=[];
intraMotor=[];
interregional=[];
intraPFCidx=[];
intraMotoridx=[];
interregionalidx=[];
intraPFCp=[];
intraMotorp=[];
interregionalp=[];
for i=1:length(Epairs)
    switch region{Epairs(i,1)}
        case 'Motor'
            switch strcmp(region(Epairs(i,1)),region(Epairs(i,2)))
                case 1
                    intraMotor=cat(2,intraMotor,normrho(:,Epairs(i,3)));
                    intraMotoridx=cat(1,intraMotoridx,Epairs(i,:));
                    intraMotorp=cat(2,intraMotorp,pval(:,Epairs(i,3)));
                    
                case 0
                    interregional=cat(2,interregional,normrho(:,Epairs(i,3)));
                    interregionalidx=cat(1,interregionalidx,Epairs(i,:));
                    interregionalp=cat(2,interregionalp,pval(:,Epairs(i,3)));
                    
            end
        case 'PFC'
            switch strcmp(region(Epairs(i,1)),region(Epairs(i,2)))
                case 1
                    intraPFC=cat(2,intraPFC,normrho(:,Epairs(i,3)));
                    intraPFCidx=cat(1,intraPFCidx,Epairs(i,:));
                    intraPFCp=cat(2,intraPFCp,pval(:,Epairs(i,3)));
                    
                case 0
                    interregional=cat(2,interregional,normrho(:,Epairs(i,3)));
                    interregionalidx=cat(1,interregionalidx,Epairs(i,:));
                    interregionalp=cat(2,interregionalp,pval(:,Epairs(i,3)));
                    
            end
    end
    
end

%% Trial z-score
% define trial
Events=CueStimulusTimes(RGC);
trial=RGC;
fq=[2:2:30 32:5:80 80:10:200];
prestim=0.2;
poststim=1.5;

trial(badRT)=[];
Events(badRT)=[];
baseline=ITIStimulusTimes(trial-1);
baseline_dur=repmat(min(ITIStimulusDuration(trial)),length(trial),1);
[zscore,pchng]=trial_zscore2(filt,fs,fq,Events,prestim,poststim,baseline,baseline_dur);
[zscoreb,pchngb]=trial_zscore2(refbi,nfs,[10 400],Events,prestim,poststim,baseline,baseline_dur);
parfor_progress(length(labels));
parfor i=1:length(labels)
[zscorec(:,:,:,i),pchngc(:,:,:,i)]=trial_zscore2(CAR(:,i),fs,fq,Events,prestim,poststim,baseline,baseline_dur);
parfor_progress;
end
parfor_progress(0);

r=ceil(size(filt,2)/16);
k=1;
for ii=1:r
figure
for i=1:16
subplot(4,4,i)
imagesc(-prestim:1/fs:poststim,fq(1:end-1),squeeze(mean(zscore(:,:,:,k),3)));set(gca,'Clim',[-1.5 1.5 ],'Ydir','normal');
title(labels(k))
k=k+1;
end
colormap jet

end

suptitle('Common Avg-silent dist')

%% ERPAC
Events=RGCR;
trial=RGC;
prestim=1;
poststim=1;
TrialType=TrialtypeR;
trial(badRT)=[];
Events(badRT)=[];
baseline=ITIStimulusTimes(trial-1);
baseline_dur=repmat(min(ITIStimulusDuration(trial)),length(trial),1);
[complex]=packagedata(hilbert(eegfilt(CAR(:,elec)',fs,4,8)'),fs,Events,prestim,poststim,baseline,baseline_dur);
complex.trials=permute(complex.trials,[2 3 1]);
complex.baseline=permute(complex.baseline,[2 3 1]);

[pac1W,pac2W]=erpacAA(complex.trials(:,find(TrialType==1),:),W,Epairs);
[pac1L,pac2L]=erpacAA(complex.trials(:,find(TrialType==2),:),L,Epairs);
[pac1S,pac2S]=erpacAA(complex.trials(:,find(TrialType==3),:),S,Epairs);