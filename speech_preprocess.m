%% import 30k data
Path=[pwd '/'];
FileName='datafilePH101816_LEFT_MER0001';
% 6 contact
ch=9:14;
%32 contact 
ch= 5:32;
%54 contact
ch=[5:32 128+(1:26)];

stripsz=length(ch);
load('/home/richardsonlab/Dropbox/Filters/30kfilters2015b.mat')
parfor_progress(length(ch));
parfor i=1:length(ch)
[~, raw, IDs] = GetAnalogData([FileName '.ns5'], 30e3,ch(i));

raw=bsxfun(@minus,raw,mean(raw,1));

tmp=unpowerline2(raw,30e3,1);

tmp=filtfilt(low30k,tmp);
tmp=downsample(tmp,5);
tmp=filtfilt(low6k,tmp);
tmp=downsample(tmp,5);
filt(:,i)=tmp;
parfor_progress;

end
parfor_progress(0);
labels=arrayfun(@(x) ['Strip' num2str(x)],1:size(filt,2),'UniformOutput',0);
addpath(genpath('/home/tom/Documents/MATLAB/eeglab13_4_4b'))
rmpath(genpath('/home/tom/Documents/MATLAB/eeglab13_4_4b/functions/octavefunc/'))

eegplot(eegfilt(filt',1200,1,0),'srate',1200)
nfs=1200;
clearvars raw tmp

%% define the bad channels
V=var(filt);
figure;
scatter(1:length(ch),V)
badch=[1 2 17 18];
%% referencing

%CAR
signal=bsxfun(@minus,filt,mean(filt(:,setdiff(1:length(ch),badch)),2));
eegplot(eegfilt(signal',1200,1,0),'srate',nfs)

%% ICA
[coeff,score,latent,tsquared,explained,~] = pca(filt(:,setdiff(1:length(ch),badch)));
plot(explained)
[weights,spheres]=runica(filt(:,setdiff(1:length(ch),badch))','extended',1,'pca',15,'maxsteps',1e3);
eegplot(weights*spheres*filt(:,setdiff(1:length(ch),badch))','srate',1200,'winlength',30,'dispchans',20)
[clean_data]=ahmad_icaproj(weights,spheres,filt(420*nfs:end,5:end)',[1 ]);
eegplot(ahmad_icaproj(weights,spheres,filt',[1 3]),'srate',1200,'winlength',10,'dispchans',20)


save('tmp.mat','-v7.3')


%% import NO data
mat_files = dir('*.mat');
input={'CRAW', 'CMacro_LFP'};
output={'Raw','LFP'};
[Rec]=ReadNeuroOmega_AA(input, output,mat_files);

%% Remove spike artifacts from micro

elec=3;
Ridx=1;
sr=Rec(Ridx).Raw.sr;
ts=round(Timestamps*sr);
ts(ts+1000>length(Rec(Ridx).Raw.ts))=[];
ts(ts-1000<1)=[];

ind=reshape(bsxfun(@plus,[-1000:1:1000]',ts'),[],1);
time=[-1000:1:1000]./sr; backtime=time(end:-1:1); AllTime=[time backtime];
figure('windowstyle','docked')
x=reshape(Rec(Ridx).Raw.ts(elec,ind),2001,[])';
r=mean(x); rPlus=r+1.645*std(x)/sqrt(size(x,1)); rMinus=r-1.645*std(x)/sqrt(size(x,1));
rm=min(rMinus); rp=max(rPlus);
backrMinus=rMinus(end:-1:1); rAllError=[rPlus backrMinus];
fill(AllTime,rAllError,[0 1 .45],'LineStyle','none'); hold on
plot(time,r,'Color',[0 0 0],'LineWidth',1); line([0 0],[rm rp],'Color','r')
ax=gca; ax.XLim=[min(time) max(time)];  ax.YLim=[rm rp];
%%
prespk=30; % samples used for the interpolation prior to spike
postspk=30;  % samples used for the interpolation after spike
spk=[40 220]; % samples being interpolated
spkt=Rec(Ridx).Raw.ts(elec,:);
for i=1:length(ts)

    spkt(ts(i)-spk(1):ts(i)+spk(2))= interp1([ts(i)-(spk(1)+prespk):ts(i)-(spk(1)+1)  ts(i)+(spk(2)+1):ts(i)+(spk(2)+postspk)],spkt([ts(i)-(spk(1)+prespk):ts(i)-(spk(1)+1)  ts(i)+(spk(2)+1):ts(i)+(spk(2)+postspk)]),ts(i)-spk(1):ts(i)+spk(2),'pchip');

end
x=reshape(spkt(ind),2001,[])';
r=mean(x); rPlus=r+1.645*std(x)/sqrt(size(x,1)); rMinus=r-1.645*std(x)/sqrt(size(x,1));
rm=min(rMinus); rp=max(rPlus);
backrMinus=rMinus(end:-1:1); rAllError=[rPlus backrMinus];
fill(AllTime,rAllError,[1 1 .45],'LineStyle','none'); hold on
plot(time,r,'Color',[0 0 0],'LineWidth',1);

%%
Rec(Ridx).Raw.cleants(elec,:)=spkt;

% if no spikes on the channel
elec=2;
Ridx=2;
Rec(Ridx).Raw.cleants(elec,:)=Rec(Ridx).Raw.ts(elec,:);
save('raw_NO_RC021616.mat','Rec','-v7.3')

%% preprocess NO data 
load('/home/richardsonlab/Dropbox/Filters/44kfilters.mat')
filt_micro=cell(size(Rec));filt_macro=filt_micro;
filt_micronspk=[];
% Raw data
parfor i =1: length(Rec)
    figure
   tmp=bsxfun(@minus,Rec(i).Raw.ts,mean(Rec(i).Raw.ts,2))';
   plot_fft(Rec(i).Raw.ts(1,:),Rec(i).Raw.sr);
   hold on
   tmp=unpowerline2(tmp,Rec(i).Raw.sr,1);
      plot_fft(tmp(1,:),Rec(i).Raw.sr);

   tmp=filtfilt(lp1200_44k,tmp);
   
   ds=downsample(tmp,11);
   tmp=filtfilt(lp500_4k,ds);
  plot_fft(tmp(1,:),Rec(i).Raw.sr/11);

   ds=downsample(tmp,2);
%    filt_micro{i}=filtfilt(hp2_2k,ds);
filt_micro{i}=ds;
end
srmic=2e3;

% Raw data
parfor i =1: length(Rec)
    figure
      tmp=bsxfun(@minus,Rec(i).Raw.cleants,mean(Rec(i).Raw.cleants,2))';
    plot_fft(Rec(i).Raw.cleants(1,:),Rec(i).Raw.sr);
    hold on
    tmp=unpowerline2(tmp,Rec(i).Raw.sr,1);
    plot_fft(tmp(1,:),Rec(i).Raw.sr);
    tmp=filtfilt(lp1200_44k,tmp);
    ds=downsample(tmp,11);
    tmp=filtfilt(lp500_4k,ds);
    plot_fft(tmp(1,:),Rec(i).Raw.sr/11);
    
    ds=downsample(tmp,2);
    %    filt_micro{i}=filtfilt(hp2_2k,ds);
    filt_micronspk{i}=ds;
end

srmic=2e3;
% Macro LFP data

for i =1: length(Rec)
   tmp=bsxfun(@minus,Rec(i).LFP.ts,mean(Rec(i).LFP.ts,2))';
   tmp=unpowerline2(tmp,Rec(i).LFP.sr,1);
%    tmp=filtfilt(lpFilt,tmp);
%    filt_macro{i}=filtfilt(hpFilt,tmp); 
   filt_macro{i}=resample(tmp,1200,1375)';
%       filt_macro{i}=tmp;

end
% srmac=Rec(i).LFP.sr;
% get force data
% for i =1: length(Rec)
% F{i}=resample(Rec(i).Force.ts(1:2,:)',8727,11);
% end
%% Import Audio

[~, Audio] = GetAnalogData([FileName '.ns5'], 30e3, [10269 ]);

%% Align ECOG and NO 

[EventTimesTrellis] = GetEventData([ FileName '.nev']); % trellis time stamps
plot(EventTimesTrellis*[1 1],ylim)
% visually identify T0 / T1 in Trellis as the start / end of segment of interest
T0 = [4000 2000 800];
for idx=1:length(Rec)
ET{idx} = sort([Rec(idx).DigUp Rec(idx).DigDown]);  % NO event times

Event0 = EventTimesTrellis(find(EventTimesTrellis>T0(idx),1,'first'));
tstart(idx)= Event0 - ET{1,idx}(1);
try
tend(idx) = tstart(idx) + length(Rec(idx).Raw.ts)/Rec(idx).Raw.sr; % 5 secs after last time stamp
E{idx}=filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);
catch
tend(idx) = tstart(idx) + ET{1,idx}(end)+10; % 5 secs after last time stamp
E{idx}=filt(round(tstart(idx)*nfs):round(tend(idx)*nfs),:);    
end
A{idx}=Audio(round(tstart(idx)*30e3):round(tend(idx)*30e3));
end
%%  save into task files
subj='BA080416';
save([subj '_raw_intraop23'],'E','filt_micro','filt_micronspk',...
    'filt_macro','FileName','nfs','srmic',...
    'ET','filt','EventTimesTrellis','Rec','labels','Audio','-v7.3')
tasks={'wordlist4_LT','wordlist3_LT','wordlist2_LT'};

for i=1:length(tasks)
    try 
    microns=filt_micronspk{i};
    catch
        microns=[];
    end
    microraw=filt_micro{i};

    macro=filt_macro{i};
    Ecog=E{i};
    EventTimes=ET{i}; 
    try
    RT=Rec(i).ResponseTimes;
    catch
        RT=[];
    end
    Audio=A{i};
    save([subj '_' tasks{i}],'microns','microraw','macro','Ecog','Audio',...
        'nfs','srmic','RT','labels',...
        'EventTimes','FileName','badch','-v7.3')
    if exist('AT')
       audiotimes=AT{i};
       save([subj '_' tasks{i}],'audiotimes','-append')
        
    end
end

