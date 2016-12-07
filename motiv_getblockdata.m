function [datastruct,rawdata] = motiv_getblockdata(datastruct)

info = datastruct.info;

% ID = [1:64,129:160,257:288];
% ch = ID(1:length(labels));
clc
% Import Filters
try
   filters = load([info.projroot,filesep,'Filters/filters.mat']); %preloaded filter set
   hp2500 = filters.hp2500;
   low30k = filters.low30k;
   lp10k  = filters.lp10k;
catch
    %make filters here: low30k, lp10k,
end
%% Force = .ns2, ECOG = .ns5
nev = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.nev')));
ns2 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns2')));
ns5 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns5')));

fsForce = 1000; %sampling rate for Force data
fsECOG  = 30000; %sampling rate for ECOG data
dbstop if error
count = 1;
for i = 1:length(info.datfiles)
    for j = 1:length(info.blocknames{i})
        beep off
        disp(['Running File ' num2str(info.datfiles(i)) ', Block ' num2str(j) '....'])
        disp(['Getting Force Data.....'])
%         try
%                 GetAnalogData(fullfile(info.trellispath,ns2{info.datfiles(i)}),...
        [~, Force, AnalogElectrodeIDs] = ...
                GetAnalogData(fullfile(info.trellispath,ns2{i}),...
                fsForce, [10241 10242], str2double(info.tstart{i}(j)),str2double(info.tstop{i}(j)));
%         catch
%         end
        % Remove 10241 and 10242 channels
        if find(AnalogElectrodeIDs==10241)>1
            AnalogElectrodeIDs = AnalogElectrodeIDs(1:find(AnalogElectrodeIDs==10241)-1);
        elseif find(AnalogElectrodeIDs==10241)==1
            AnalogElectrodeIDs = AnalogElectrodeIDs(3:end);
        end
        disp(['Running File ' num2str(info.datfiles(i)) ', Block ' num2str(j) '....'])
        disp(['Getting ECOG Data.....'])

% addAttachedFiles(pool,{'low30k','lp10k','hp2500'})
tic
parfor k=1:length(AnalogElectrodeIDs)
% for k=1:length(AnalogElectrodeIDs)
    disp(['Loop ' num2str(k) ' of ' num2str(length(AnalogElectrodeIDs))])
    [~, rawdata, ~] = ...
        GetAnalogData(fullfile(info.trellispath,ns5{i}),...
        fsECOG, AnalogElectrodeIDs(k), str2double(info.tstart{i}(j)),str2double(info.tstop{i}(j)));
%     GetAnalogData(fullfile(info.trellispath,ns5{info.datfiles(i)}),...
    rawdata=bsxfun(@minus,rawdata,mean(rawdata,1));
    
    rawdata=unpowerline(rawdata,fsECOG,2);
    f1=filtfilt(low30k,rawdata);
    
    ds1=downsample(f1,3);
    
    f2=filtfilt(lp10k,ds1);
    ds2=downsample(f2,4);
    filtered(:,k)=filtfilt(hp2500,ds2);
    
end
toc
beep on
beep
datastruct.filtdata{count}  = filtered;
datastruct.forcedata{count} = Force;

clear filtered f1 f2 rawdata Force ds1 ds2
count = count+1;
disp(['Count is ' num2str(count)])

% filtdata = datastruct.filtdata{count-1};
% forcedata = datastruct.forcedata{count-1};
% save(['File' num2str(info.datfiles(i)) '_Block' char(info.blocknames{i}(j-1)) '_data.mat'],'filtdata','forcedata','-v7.3')
% clear filtdata forcedata filtered
% save('datastruct.mat','datastruct','-v7.3')
    end
end
datastruct.AnalogElectrodeIDs   = AnalogElectrodeIDs;
datastruct.fsECOGfilt           = fsECOG/12;
datastruct.fsForce              = fsForce;

% count=2;
% eegplot(datastruct.filtdata{2}','srate',fsECOG/12)
% eegplot(datastruct.filtdata{2}','srate',datastruct.fsECOGfilt)
% clearvars f1 f2 dsdata rawdata
% nfs=fsECOG/12;
% save(fullfile(info.patientdir,'raw.mat'),'raw','-v7.3')

% if ~exist(fullfile(info.patientdir,'datastruct.mat'),'file')
    save(fullfile(info.patientdir,'datastruct.mat'),'datastruct','-v7.3')
% else
%     for i = 1:10
%         if ~exist(fullfile(info.patientdir,['datastruct0' num2str(i) '.mat']),'file');
%             save(fullfile(info.patientdir,['datastruct0' num2str(i) '.mat']),'datastruct','-v7.3');
%         end
%     end
% end
% end
% end
end
%
% %Load data
%         try
%         [tForce, Force, AnalogElectrodeIDs] = GetAnalogData(fullfile(Trellispath,ns2{info.datfiles(i)}), fsForce, [10241 10242],str2double(tstart{j}),str2num(tstop{j}));
%         catch
%         [tForce, Force, AnalogElectrodeIDs] = GetAnalogData(fullfile(Trellispath,ns2{i}), fsForce, [10241 10242],str2double(tstart{j}),str2num(tstop{j}));
%         end
%         %[tLFP, LFP, ~, lfpFileInfo] = GetAnalogData_sk([LFPfname,LFPfext], fs, ch);
%         %  ch=AnalogElectrodeIDs(1:end-2);
%         if AnalogElectrodeIDs(end)==10242
%             ID = AnalogElectrodeIDs(1:end-2);
%         else
%             ID = AnalogElectrodeIDs(3:end);
%         end
%         try
%         [tECOG, ECOGraw, ID] = GetAnalogData(fullfile(Trellispath,ns2{info.datfiles(i)}), fsECOG, ID,str2num(tstart{j}),str2num(tstop{j}));
%         catch
%         [tECOG, ECOGraw, ID] = GetAnalogData(fullfile(Trellispath,ns2{i}), fsECOG, ID,str2num(tstart{j}),str2num(tstop{j}));
%         end
%         %%
%         ECOGraw = bsxfun(@minus,ECOGraw,mean(ECOGraw,1));
%         disp('unpowerline....')
%         ECOGraw = unpowerline2(ECOGraw,fsECOG,2);
%
% %         disp('filter data...')
% % ECOGds    = filtfilt(low30k,ECOGraw);
% % ECOGds = downsample(ECOG,3); %10000
% % ECOGds = filtfilt(lp10k,ds);
% % ECOGds = downsample(ECOG,4);
% % ECOGds = filtfilt(hp2500,ds);
%
% % ECOGdsrate = ECOGrawsrate/12;
% % eegplot(ECOG','srate',fsECOG)
% eval(['datastruct.tstart.block' num2str(block(i)) ' = tstart;']);
% eval(['datastruct.tstop.block' num2str(block(i)) ' = tstop;']);
%
% try
% eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.IDs  = [10241 10242];']);
% eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.Force= Force;']);
% eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.srate= fsForce;']);
% eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.time = tForce;']);
% catch
% eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) ' = [];']);
% end
%
% try
% eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.IDs       = ID;']);
% eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.rawdata   = ECOGraw;']);
% eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.rawsrate  = fsECOG;']);
% % eval(['datastruct.1kdata.block' num2str(info.reply(block(i))) '.session' num2str(j) '.ECOGds    = ECOGds;']);
% % eval(['datastruct.1kdata.block' num2str(info.reply(block(i))) '.session' num2str(j) '.ECOGdsrate= fsECOG;']);
% eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.ECOGtime  = tECOG;']);
% catch
% eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) ' = [];']);
% end
%     end
%     end
%
% end
%
