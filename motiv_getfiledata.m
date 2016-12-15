function [datastruct] = motiv_getfiledata(info)

if ~isfield(info,'nev')
    info.nev = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.nev')));
    info.ns2 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns2')));
    info.ns5 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns5')));
else
end
datastruct.ptid = info.patientname;
datastruct.info = info;
datastruct.datfiles = info.ns5;
% ID = [1:64,129:160,257:288];
% try
%     ch = ID(1:length(info.labels));
% catch
%     ch = ID(:);
% end
clc
% Import Filters
try
   filters = load([info.projroot,filesep,'Filters/30kfilters2015b.mat']); %preloaded filter set
   hp1200 = filters.hp1200;
   low30k = filters.low30k;
   low6k  = filters.low6k;
catch
    %make filters here: low30k, lp10k,
end

fsForce = 1000; %sampling rate for Force data
fsECOG  = 30000; %sampling rate for ECOG data
dbstop if error
count = 1;
for i = 1:length(info.ns5)
    beep off
    disp(['Running File ' num2str(info.ns5{i}) '....'])
    disp('Getting Force Data.....')
    [~,Force,AnalogElectrodeIDs] = GetAnalogData(...
        fullfile(info.trellispath,info.ns2{i}),fsForce,[10241 10242]);
    % Remove 10241 and 10242 channels
    if find(AnalogElectrodeIDs==10241)>1
        AnalogElectrodeIDs = AnalogElectrodeIDs(1:find(AnalogElectrodeIDs==10241)-1);
    elseif find(AnalogElectrodeIDs==10241)==1
        AnalogElectrodeIDs = AnalogElectrodeIDs(3:end);
    end
    disp(['Running File ' num2str(info.ns5{i}) '....'])
    disp(['Getting ECOG Data.....'])
    
tic
parfor_progress(length(AnalogElectrodeIDs));
parfor k=1:length(AnalogElectrodeIDs)
    disp(['Loop ' num2str(k) ' of ' num2str(length(AnalogElectrodeIDs))])
    [~, rawdata, ~] = GetAnalogData(...
        fullfile(info.trellispath,info.ns5{i}),fsECOG, AnalogElectrodeIDs(k));

    rawdata=bsxfun(@minus,rawdata,mean(rawdata,1));
    rawdata=unpowerline(rawdata,fsECOG,2);
    f1=filtfilt(low30k,rawdata);
    
    ds1=downsample(f1,5);
    f2=filtfilt(low6k,ds1);
    ds2=downsample(f2,5);
    filtered(:,k)=filtfilt(hp1200,ds2);
    
end
toc
parfor_progress(0);
beep on
beep
datastruct.forcedata{count} = Force;
datastruct.filtdata{count}  = filtered;

clear filtered f1 f2 rawdata Force ds1 ds2
count = count+1;
disp(['Count is ' num2str(count)])

end
datastruct.IDs   = AnalogElectrodeIDs;
end

% strfind(datastruct.info.ns5,'_')
% datastruct.fileNames{1}=datastruct.info.ns5{1}(17+1:end-4)
% datastruct.fileNames{2}=datastruct.info.ns5{2}(17+1:end-4)