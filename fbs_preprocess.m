%% Select Channels
ch=[[1:4],[9:16]]; choi = ch;
%ch=1:32;

%% Initialize Force and LFP files
clear pID
pID = cell2mat(inputdlg('What is patient ID?','pID',1,{' '}));
%---------------------------------------------------------------------------------------------------------------
[Ffname, Fpathname, ~] = uigetfile({'.ns2'}, 'Choose Force File'); 
[Fpathname,Ffname,Ffext] = fileparts(fullfile(Fpathname,Ffname));
[LFPfname, LFPpathname, ~] = uigetfile({'*.ns2','ns2';'*.ns5','ns5'}, 'Choose LFP File');
[LFPpathname,LFPfname,LFPfext] = fileparts(fullfile(LFPpathname,LFPfname));

%Warn if fnames don't match
if strcmp(Ffname,LFPfname)~=1
    choice = questdlg(sprintf('Did you choose the same Force and LFP files??\n Force - %s\n LFP - %s',Ffname,LFPfname), ...
	'Yes', 'No')
    switch choice
        case 'No'
            return
    end
end

%Set fs
if strcmp(Ffext,'.ns2')==1
    fsF=1000;
else
    warndlg('Force file should be .ns2!!')
    return
end
if strcmp(LFPfext,'.ns2')==1
    fs=1000;
elseif strcmp(LFPfext,'.ns5')==1
    fs=30000;
end

%Load data
[tF, Force, AnalogElectrodeIDs] = GetAnalogData(fullfile(Fpathname,[Ffname,Ffext]), fsF, [10241 10242]);
%[tLFP, LFP, ~, lfpFileInfo] = GetAnalogData_sk([LFPfname,LFPfext], fs, ch);
 ch=AnalogElectrodeIDs(1:end-2);
[tLFP, LFP, ~] = GetAnalogData(fullfile(LFPpathname,[LFPfname,LFPfext]), fs, ch);

Force=single(Force);
LFP=single(LFP);

%Force = Force(:,2:-1:1);  %Reverse handles for MM

%% Cut Recording
tF1 = tF./(60);
f2=figure(2); clf, plot(tF1, Force);

tstart = 7*60;  tend = 14.1*60;

t_thr = round([tstart, tend] .*fs );

Force = Force(t_thr(1):t_thr(2),:);
LFP = LFP(t_thr(1):t_thr(2),:);
tF = tF(t_thr(1):t_thr(2)) - t_thr(1)/fs;

f2=figure(2); clf, plot(tF, Force);

%% Downsample if .ns5
if fs==30000       %Downsample
[LFP, fs] = stath_decimate (LFP, fs, 10);
[LFP, fs] = stath_decimate (LFP, fs, 3);
tLFP=linspace(tLFP(1),tLFP(end),size(LFP,1));
end

%% Response Times
Lskip=0; Rskip=0;
[~, LR, ~] = MotivationTaskResponseTrig(Force(:,1), 100, 0.5, 0.5, 1000, 300,1);
[~, RR, ~] = MotivationTaskResponseTrig(Force(:,2), 100, 0.5, 0.5, 1000, 300,1);
LR=LR/fsF; LR(1:Lskip)=[]; RR=RR/fsF; RR(1:Rskip)=[];%Modulate based on skips....
R=[RR;LR];
clc

f2=figure(2); clf, plot(tF, Force);
hold on
line([LR,LR]',repmat([-50,100],length(LR),1)','linestyle','--','color','b')
line([RR,RR]',repmat([-50,100],length(RR),1)','linestyle','--','color','r')
%% Text Data
clear('Bilateral*'), pause(0.1)
uiimport('-file')

a=who('Bilateral*');
eval(['txtdata=' a{1} ';'])

%% Events
[EventTimes] = GetEventData(fullfile(LFPpathname,[LFPfname '.nev']));
if exist('t_thr','var'), EventTimes=EventTimes - t_thr(1)/fs; EventTimes(EventTimes<0)=[]; end

f2=figure(2); clf, plot(tF, Force); yl=ylim;
E=EventTimes; E(1:3)=[]; E=downsample(E,4);
hold on; plot(EventTimes*[1 1], [yl(1), -25], '--k'); plot(E*[1,1], [0, yl(2)], '--r')
for i=1:length(EventTimes)
   text(EventTimes(i),-20,sprintf('%d',i))
end
uiwait(f2)

%Skips
prompt = {'Event skips:','Trial Skips:'};
dlg_title = 'Skips';
num_lines = 1;
def = {'1:3','0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
eskip=answer{1};
tskip=str2double(answer{2});

eval(['EventTimes([' eskip '])=[];'])

if tskip>0
    txtdata(1:tskip,:)=[];
end

%% Trial
[CueStimulusTimes, CueStimulusDuration, CueTrial, CueTrialNum] = MotivationTaskStimulusBlock2(EventTimes, 0, 4, 1, txtdata);
[CommandStimulusTimes, CommandStimulusDuration, CommandTrial, CommandTrialNum] = MotivationTaskStimulusBlock2(EventTimes, 0, 4, 2, txtdata);
[FeedbackStimulusTimes, FeedbackStimulusDuration, FeedbackTrial, FeedbackTrialNum] = MotivationTaskStimulusBlock2(EventTimes, 0, 4, 3, txtdata);
[ITIStimulusTimes, ITIStimulusDuration, ITITrial, ITITrialNum] = MotivationTaskStimulusBlock2(EventTimes, 0, 4, 4, txtdata);
[LRStimulusTimes, LRStimulusDuration, LRTrial, LRTrialNum] = MotivationTaskStimulusBlock2(LR, 0, 1, 1, txtdata);
[RRStimulusTimes, RRStimulusDuration, RRTrial, RRTrialNum] = MotivationTaskStimulusBlock2(RR, 0, 1, 1, txtdata);
[RStimulusTimes, RStimulusDuration, RTrial, RTrialNum] = MotivationTaskStimulusBlock2(R, 0, 1, 1, txtdata);

LRTrialNum=ITITrialNum; RRTrialNum=ITITrialNum; %Correct for error in LFPStimulus....
trial=CueTrial;
trial.RR=RR; trial.LR=LR; trial.R=R;

%% Labels

%Labels for DBS
labels={};
for i=1:size(LFP,2)
   labels{i,1}=sprintf('C%.2d',i); 
end

%% Save Trial
svname=[pID '_trial'];
%save(svname,'trial','Force','CueStimulusTimes','CommandStimulusTimes','FeedbackStimulusTimes','ITIStimulusTimes','RRStimulusTimes','LRStimulusTimes','RStimulusTimes','CueTrialNum','CommandTrialNum','FeedbackTrialNum','ITITrialNum','RTrialNum','RRTrialNum','LRTrialNum')
save(svname,'trial','Force','CueStimulusTimes','CommandStimulusTimes','FeedbackStimulusTimes','ITIStimulusTimes')

%% Save LFP
raw=single(LFP);
raw=raw(:,1:length(labels));

svnameLFP=[pID '_raw'];
save(svnameLFP,'raw','fs','labels','-v7.3') 

%% Select channels and Apply Unpowerline
%choi = ch;
referenced = unpowerline(raw(:,choi),fs,.5);

svnameLFP = [pID '_noref'];
save(svnameLFP,'referenced','fs','labels','-v7.3')

%% Modify Labels

labels={};
for i=1:4
    labels{i}=sprintf('STN_%d',i);
end
for i=5:12
    labels{i}=sprintf('Cortex_%d',i-4);
end
labels = labels';
save(svnameLFP,'-append','labels')

%% Segment

s=1:size(referenced,2);
LFP=referenced(:,s); clear referenced
labels=labels(s);

trialdata={};
trialforce={};


%-------------------% Set Stimulus Times------------------------
centers = {'Cue','Command','Response','Cue'};
tcntr0 = [0,2;0,2;-1,1];
tcntr = {[0*fs+1:2*fs], [0*fs+1:2*fs], [-1*fs+1:1*fs], [-1*fs+1:0*fs]};
for cntr_i=1:length(centers)
    cntr=centers{cntr_i};
    pre=3; preT=pre;
    post=3; postT=post; 
    RTthr=max(ITIStimulusTimes-CommandStimulusTimes); %Max time we will look for a response is until Feedback happens
    
    if strcmp(cntr,'Response')~=1
        eval(['stim=' cntr 'StimulusTimes;'])
    elseif strcmp(cntr,'Response')==1
        stim=nan(1,length(CommandStimulusTimes));
        for i=1:length(CommandStimulusTimes)
            if ismember(i,trial.Right)
                tdiff=trial.RR-CommandStimulusTimes(i);
                tdiff=min(intersect(tdiff(tdiff>0),tdiff(tdiff<RTthr)));
                if ~isempty(tdiff)
                    stim(i)=CommandStimulusTimes(i)+tdiff;
                else
                    stim(i)=nan;
                end
            elseif ismember(i,trial.Left)
                tdiff=trial.LR-CommandStimulusTimes(i);
                tdiff=min(intersect(tdiff(tdiff>0),tdiff(tdiff<RTthr)));
                if ~isempty(tdiff)
                    stim(i)=CommandStimulusTimes(i)+tdiff;
                else
                    stim(i)=nan;
                end
            end
        end
    end
    
    %-------------------------Set Trial Data - ----------------------------
    %tt=(round(-pre*fs)+1 : round(post*fs))/fs;
    tt=tcntr{cntr_i};
    trialdata{cntr_i}=nan(length(tt),length(stim),length(labels),'single');
    trialforce{cntr_i}=nan(length(tt),length(stim),2,'single');
    
    for ch=1:size(labels,1)
        for s=1:size(stim,2)
            if ~isnan(stim(s))
                stimdata                    = LFP(round(tt+stim(s)*fs),ch);
                trialdata{cntr_i}(:,s,ch)   = stimdata;
                stimforce                   = Force(round(tt+stim(s)*fs),:);
                trialforce{cntr_i}(:,s,:)   = stimforce;

            end
        end
    end
    
end

save(svname,'-append','trialdata','trialforce','labels')
