%%                   Initialize Variables
%%%%%%%%%%%%%%%%
% addpath('/Users/markrichardson/Documents/Projects/Motivation')
clear
clc
info = motiv_prepinfo('RH021016');
info = motiv_getblockinfo(info);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datastruct.PtId = info.patientdir(max(strfind(info.patientdir,'/'))+1:end);
datastruct.info = info;
[datastruct] = motiv_getblockdata(datastruct); % block=3; ID=[]; 

%%
block=2;
sesh=1;
data = eval(['datastruct.k1data.block' num2str(block) '.session' num2str(sesh) '.rawdata']);
figure;
for ch=1:size(data,2)
hold on; plot(data(:,ch)-2500*ch)
end
% for ch=2:32
% hold on; plot(data(:,ch)-data(:,ch-1)-2500*ch)
% end
%%
%USER INPUT: Choose Block.................
% for i = 1:length(info.reply)
block = 3; 
ID    = []; %Grapevine, see GetAnalogData
load '/Users/markrichardson/Documents/Projects/Motivation/Filters/filters.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ECOGraw = eval(['datastruct.block' num2str(datastruct.info.reply(block) '.ECOGraw'])

raw_srate       = 30000;
ECOGraw         = filtfilt(low30k,ECOGraw); %low_pass
ECOGrs          = resample(ECOG,3); %10000
ECOG = filtfilt(lp10k,ds);
ds = downsample(ECOG,4);
ECOG = filtfilt(hp2500,ds);

fsECOG = fsECOG/12;
eegplot(ECOG','srate',fsECOG)

eval(['datastruct.block' num2str(info.reply(block)) '.IDs       = ID;']);
eval(['datastruct.block' num2str(info.reply(block)) '.ECOG      = ECOG;']);
eval(['datastruct.block' num2str(info.reply(block)) '.ECOGsrate = fsECOG;']);
eval(['datastruct.block' num2str(info.reply(block)) '.ECOGtime  = tECOG;']);
eval(['datastruct.block' num2str(info.reply(block)) '.Force     = Force;']);
eval(['datastruct.block' num2str(info.reply(block)) '.Forcesrate= fsForce;']);
eval(['datastruct.block' num2str(info.reply(block)) '.Forcetime = tForce;']);

keep datastruct 

save(fullfile(datastruct.info.projroot,'datastruct'),'datastruct','-v7.3')



%% BIPOLAR REFERENCING
% All depth macroelectrodes (depth) should be bipolared. Micro electrodes should not be
% bipolared. 
m=8; % 8/8 v. 10/6 where micro = 6;
% elec={'RAMY', 'RHIPP','RPT', 'LAMY', 'LHIPP', 'LBT', 'LPSTG', 'LSTS', 'LA','LSM','LSI','LMI','LFO','LOF','LMF' };
signal=zeros(size(ECOG));
temp2=[];
micro=[];
for e=1:length(elec)
    temp=find(~cellfun(@isempty,strfind(labels,elec{e})));
    if length(temp)==16
         micro=cat(1,micro,temp(m+1:end));
        temp=temp(1:16-m);
    end
    signal(:,temp(1:end-1))=bsxfun(@minus,ECOG(:,temp(1:end-1)),ECOG(:,temp(2:end)));
    temp2=cat(1,temp2,temp);
end
signal(:,micro)=ECOG(:,micro);

% Behnke Fried 
temp=find(~cellfun(@isempty,strfind(labels,'LBF')));
signal(:,temp)=bsxfun(@minus,ECOG(:,temp),ECOG(:,temp(length(temp))));
temp2=cat(1,temp2,temp);

%% COMMON AVERAGE REFERENCING

temp=setdiff(1:size(ECOG,2),temp2);
signal(:,temp)=bsxfun(@minus,ECOG(:,temp),mean(ECOG(:,temp),2));


eegplot(signal','srate',srate/12)
%% define the bad channels
if ~exist(badch)
    badch=[];
else
    ECOG(:,badch)=0;
end
%% ICA

[coeff,score,latent,tsquared,explained,~] = pca(ECOG(420*nfs:end,5:end));
plot(explained)
[weights,spheres]=runica(ECOG(420*nfs:end,5:end)','extended',1,'pca',5,'maxsteps',1e3);
eegplot(weights*spheres*ECOG(420*nfs:end,5:end)','srate',1200,'winlength',30,'dispchans',20)
[clean_data]=ahmad_icaproj(weights,spheres,ECOG(420*nfs:end,5:end)',[1 ]);
eegplot(ahmad_icaproj(weights,spheres,ECOG',[1 3]),'srate',1200,'winlength',10,'dispchans',20)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Start Psychometrics Processing                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Psychometrics
% Get Force
% Left = 10241(blue); Right = 10242(red)
srate=1000;
[time1kHz, Force, IDs] = GetAnalogData(fullfile(path2data,[FileName,'.ns2']),srate,[10241 10242]);

figure('windowstyle','docked'); plot(time1kHz, Force) % force(y) v. time
%% Two Trials
tstart= 1;
tend= round(580*nfs);
signal=signal(tstart:tend,:);
Force=Force(tstart:580*srate,:);
time1kHz=1/srate:1/srate:length(Force)/srate;

%%
% tstart2= round(650*nfs);
% signal3 = signal3(tstart2:end,:);
% [time1kHz2, Force2, IDs] = GetAnalogData([Path FileName '.ns2'], fs, [10241 10242]);
% Force2 = Force2(650*fs:end,:);


%%
% Get Left and Right Reponse Times

[LeftResponseTimes,movementLt]  = MotivationTaskResponseTrig(Force(:,1), 100, srate, 250, 100, 0.07);
[RightResponseTimes,movementRt]  = MotivationTaskResponseTrig(Force(:,2), 100, srate, 250, 100, 0.07);


figure('windowstyle','docked'); plot(time1kHz, Force)
hold on; plot((1/srate)*LeftResponseTimes*[1 1], ylim, 'b');
hold on; plot((1/srate)*RightResponseTimes*[1 1], ylim, 'g');
hold on; plot((1/srate)*movementLt.deflection'*[1 1], ylim, 'r');
hold on; plot((1/srate)*movementRt.deflection'*[1 1], ylim, 'y');
scatter(movementLt.locs/srate,movementLt.peaks,'b','filled')
scatter(movementRt.locs/srate,movementRt.peaks,'r','d','filled')


%Timestamps for task parameters -- Just gotta know the order they appear
[EventTimes] = GetEventData(fullfile(path2data,[FileName,'.nev']));
EventTimes = EventTimes(1:391,:);
% 
% [EventTimes2] = GetEventData(fullfile(Path,[FileName,'.nev']));
% EventTimes2 = EventTimes2(392:end,:);
% EventTimes2 = bsxfun(@minus,EventTimes2, 650);



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
%For events 1:4
[CueStimulusTimes, CueStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 1, celldata); %Event1
[CommandStimulusTimes,~, CommandTrial, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 2, celldata);     %Event2
[FeedbackStimulusTimes, ~, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 3, celldata);              %Event3
[ITIStimulusTimes, ITIStimulusDuration, ~, ~] = MotivationTaskStimulusBlock(EventTimes, skips, 4, 4, celldata); %Event4
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
        RGCR(i,1)=RightResponseTimes(intersect(find(RightResponseTimes> srate*CommandStimulusTimes(RGC(i))),find( RightResponseTimes < srate*CommandStimulusTimes(RGC(i))+ 2*srate)));
    catch
        RGCR(i)=0;
       
    end
end
RGC=RGC(RGCR>0);
RGCR=RGCR(RGCR>0)/srate;
rind=cell2mat(arrayfun(@(x) intersect(find(movementRt.locs> srate*x),find( movementRt.locs < srate*x+ 2*srate)),RGCR,'UniformOutput',false));
movementRt.locs=movementRt.locs(rind);
movementRt.peaks=movementRt.peaks(rind);
movementRt.deflection=movementRt.deflection(rind);

for i=1:length(LGC)
    try
        LGCR(i,1)=     LeftResponseTimes(intersect(find(LeftResponseTimes> srate*CommandStimulusTimes(LGC(i))),find( LeftResponseTimes < srate*CommandStimulusTimes(LGC(i))+ 2*srate)));
    catch
        LGCR(i)=0;
  
    end
end
LGC=LGC(LGCR>0);
LGCR=LGCR(LGCR>0)/srate;
lind=cell2mat(arrayfun(@(x) intersect(find(movementLt.locs> srate*x),find( movementLt.locs < srate*x+ 2*srate)),LGCR,'UniformOutput',false));

movementLt.locs=movementLt.locs(lind);
movementLt.peaks=movementLt.peaks(lind);
movementLt.deflection=movementLt.deflection(lind);




%%
% 
% 
% [CueStimulusTimes2, CueStimulusDuration2, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 1, celldata2);
% [CommandStimulusTimes2,~, CommandTrial2, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 2, celldata2);
% [FeedbackStimulusTimes2, ~, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 3, celldata2);
% [ITIStimulusTimes2, ITIStimulusDuration2, ~, ~] = MotivationTaskStimulusBlock(EventTimes2, skips, 4, 4, celldata2);
% switch task
%     case 1
%         RGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.Go));
%         LGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.Go));
%         RNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.NoGo));
%         LNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.NoGo));
%     case 2
%         RGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.Go));
%         RGC2=cat(2,RGC2, intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.GoFast)));
%         LGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.Go));
%         LGC2=cat(2,LGC2, intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.GoFast)));
%         RNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Right,CommandTrial2.NoGo));
%         LNGC2=intersect(CommandTrial2.Correct,intersect(CommandTrial2.Left,CommandTrial2.NoGo));
% end
% RNGC2=setdiff(RNGC2,1);
% LNGC2=setdiff(LNGC2,1);
% RGC2=setdiff(RGC2,1);
% LGC2=setdiff(LGC2,1);
% 
% for i=1:length(RGC2)
%     try
%         RGCR2(i,1)=RightResponseTimes2(intersect(find(RightResponseTimes2> fs*CommandStimulusTimes2(RGC2(i))),find( RightResponseTimes2 < fs*CommandStimulusTimes2(RGC2(i))+ 2*fs)));
%     catch
%         RGCR2(i)=0;
%        
%     end
% end
% RGC2=RGC2(RGCR2>0);
% RGCR2=RGCR2(RGCR2>0)/fs;
% rind2=cell2mat(arrayfun(@(x) intersect(find(movementRt2.locs> fs*x),find( movementRt2.locs < fs*x+ 2*fs)),RGCR2,'UniformOutput',false));
% movementRt2.locs=movementRt2.locs(rind2);
% movementRt2.peaks=movementRt2.peaks(rind2);
% movementRt2.deflection=movementRt2.deflection(rind2);
% 
% for i=1:length(LGC2)
%     try
%         LGCR2(i,1)=     LeftResponseTimes2(intersect(find(LeftResponseTimes2> fs*CommandStimulusTimes2(LGC2(i))),find( LeftResponseTimes2 < fs*CommandStimulusTimes2(LGC2(i))+ 2*fs)));
%     catch
%         LGCR2(i)=0;
%   
%     end
% end
% LGC2=LGC2(LGCR2>0);
% LGCR2=LGCR2(LGCR2>0)/fs;
% lind2=cell2mat(arrayfun(@(x) intersect(find(movementLt2.locs> fs*x),find( movementLt2.locs < fs*x+ 2*fs)),LGCR2,'UniformOutput',false));
% 
% movementLt2.locs=movementLt2.locs(lind2);
% movementLt2.peaks=movementLt2.peaks(lind2);
% movementLt2.deflection=movementLt2.deflection(lind2);
%% save
% mkdir('/home/richardsonlab/Dropbox/STN_cortical connectivity/Epilepsy/AW093014')
% save('/home/richardsonlab/Dropbox/STN_cortical connectivity/JM071415/JM08_RT.mat','-v7.3')

%%  Trial processing
for e=1:length(LGC)
%     TrialtypeL(e)=sum(ismember(CommandTrial.Win,LGC(e)))+2*sum(ismember(CommandTrial.Lose,LGC(e)))+3*sum(ismember(CommandTrial.Squeeze,LGC(e)));
    EEG.event(1).type='Cue';
    EEG.event(1).latency=srate*(CueStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    EEG.event(2).type='Command';
    EEG.event(2).latency=srate*(CommandStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    EEG.event(3).type='Squeeze';
    EEG.event(3).latency=srate*(LGCR(e)-ITIStimulusTimes(LGC(e)-1));
    EEG.event(4).type='Feedback';
    EEG.event(4).latency=srate*(FeedbackStimulusTimes(LGC(e))-ITIStimulusTimes(LGC(e)-1));
    ELeft{e}=EEG;
end
[trialsL,badL]=check_trials(signal,ITIStimulusTimes(LGC-1),ITIStimulusTimes(LGC),ELeft,nfs);

for e=1:length(RGC)
    TrialtypeR(e)=sum(ismember(CommandTrial.Win,RGC(e)))+2*sum(ismember(CommandTrial.Lose,RGC(e)))+3*sum(ismember(CommandTrial.Squeeze,RGC(e)));
    EEG.event(1).type='Cue';
    EEG.event(1).latency=srate*(CueStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    EEG.event(2).type='Command';
    EEG.event(2).latency=srate*(CommandStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    EEG.event(3).type='Squeeze';
    EEG.event(3).latency=srate*(RGCR(e)-ITIStimulusTimes(RGC(e)-1));
    EEG.event(4).type='Feedback';
    EEG.event(4).latency=srate*(FeedbackStimulusTimes(RGC(e))-ITIStimulusTimes(RGC(e)-1));
    [pk loc]=max(Force(round(RGCR(e)*srate):round(RGCR(e)*srate)+srate,2));
    acc=(pk-Force(round(RGCR(e)*srate),2))/((loc-1)/srate);
    rxt=RGCR(e)-CommandStimulusTimes(RGC(e));
    movement(e,:)=[rxt pk acc];
    ERight{e}=EEG;
    
end
[trialsR,badR]=check_trials(signal,ITIStimulusTimes(RGC-1),ITIStimulusTimes(RGC),ERight,nfs);



%% Downsample data
dsSignal=downsample(filtfilt(lpFilt,signal),25); %2500 to 250
dsrate = nfs/25;
save('Data_small.mat','dsSignal','dsrate','labels','elec',...
    'BilateralDelayMotivationTaskIntraop','EEG','EventTimes',...
    'FeedbackStimulusTimes','Force','ITIStimulusDuration',...
    'ITIStimulusDuration','CommandStimulusTimes','CommandTrial',...
    'CueStimulusDuration','CueStimulusTimes','ELeft','ERight',...
    'movement','movementLt','movementRt','nfs','fs','rind','RNGC',...
    'LGC','LGCR','lind','LNGC','loc','micro','PtId','RGC','RGCR',...
    'RightResponseTimes','skips','time1kHz','trialskips','trialsL',...
    'trialsR','TrialtypeR')

%%
