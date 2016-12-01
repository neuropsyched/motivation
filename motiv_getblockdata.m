function [datastruct] = motiv_getblockdata(datastruct)

info = datastruct.info;
for i=1:length(info.reply)
    block(i) = str2double(info.reply(i));
end
clear i
% Import Filters
try
    load([info.projroot,filesep,'Filters/filters.mat']) %preloaded filter set
catch
    %make filters here: low30k, lp10k, 
end
%% Force = .ns2, ECOG = .ns5
nev = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.nev')));
ns2 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns2')));
ns5 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns5')));
Trellispath = fullfile(info.patientdir,'Trellis');

fsForce = 1000; %sampling rate for Force data
fsECOG  = 1000; %sampling rate for ECOG data

%Load data

for i = 1:length(block)
    if ~(length(nev)<i)
    [EventTimes] = GetEventData(fullfile(info.patientdir,'Trellis',nev{i}));
    % Plot channels
    if ~isempty(EventTimes) | ~EventTimes==0
    h = figure; plot(EventTimes*[1 1], ylim, 'k')
    title(['block' num2str(block(i))])
    waitfor(h) % pause while figure is open
    disp(['Choose start and stop times'])
    tstart = strsplit(cell2mat(inputdlg('Enter start time')),','); % :)
    tstop = strsplit(cell2mat(inputdlg('Enter stop time')),',');
    % str2num(char(cellfun(@(x) (x(1:strfind(tstart{1},',')-1)), tstart,'Uniformoutput',false)))
    else
        tstart = {'0'};
        tstop = {''};
    end
    else
        tstart = {'0'};
        tstop = {''};
    end
    if ~isequal(length(tstart),length(tstop))
        error('tstart and tstop must be the same length')
    end
    for j = 1:length(tstart)
        disp(['Running Block ' num2str(block(i)) ', Session ' num2str(j) '....'])
        try
        [tForce, Force, AnalogElectrodeIDs] = GetAnalogData(fullfile(Trellispath,ns2{block(i)}), fsForce, [10241 10242],str2double(tstart{j}),str2num(tstop{j}));
        catch
        [tForce, Force, AnalogElectrodeIDs] = GetAnalogData(fullfile(Trellispath,ns2{i}), fsForce, [10241 10242],str2double(tstart{j}),str2num(tstop{j}));
        end    
        %[tLFP, LFP, ~, lfpFileInfo] = GetAnalogData_sk([LFPfname,LFPfext], fs, ch);
        %  ch=AnalogElectrodeIDs(1:end-2);
        if AnalogElectrodeIDs(end)==10242
            ID = AnalogElectrodeIDs(1:end-2);
        else
            ID = AnalogElectrodeIDs(3:end);
        end
        try
        [tECOG, ECOGraw, ID] = GetAnalogData(fullfile(Trellispath,ns2{block(i)}), fsECOG, ID,str2num(tstart{j}),str2num(tstop{j}));
        catch
        [tECOG, ECOGraw, ID] = GetAnalogData(fullfile(Trellispath,ns2{i}), fsECOG, ID,str2num(tstart{j}),str2num(tstop{j}));
        end   
        %%
        ECOGraw = bsxfun(@minus,ECOGraw,mean(ECOGraw,1));
        disp('unpowerline....')
        ECOGraw = unpowerline2(ECOGraw,fsECOG,2);
        
%         disp('filter data...')
% ECOGds    = filtfilt(low30k,ECOGraw);
% ECOGds = downsample(ECOG,3); %10000
% ECOGds = filtfilt(lp10k,ds);
% ECOGds = downsample(ECOG,4);
% ECOGds = filtfilt(hp2500,ds);

% ECOGdsrate = ECOGrawsrate/12;
% eegplot(ECOG','srate',fsECOG)
eval(['datastruct.tstart.block' num2str(block(i)) ' = tstart;']);
eval(['datastruct.tstop.block' num2str(block(i)) ' = tstop;']);

try
eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.IDs  = [10241 10242];']);
eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.Force= Force;']);
eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.srate= fsForce;']);
eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) '.time = tForce;']);
catch
eval(['datastruct.Forcedata.block' num2str(block(i)) '.session' num2str(j) ' = [];']);
end

try
eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.IDs       = ID;']);
eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.rawdata   = ECOGraw;']);
eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.rawsrate  = fsECOG;']);
% eval(['datastruct.1kdata.block' num2str(info.reply(block(i))) '.session' num2str(j) '.ECOGds    = ECOGds;']);
% eval(['datastruct.1kdata.block' num2str(info.reply(block(i))) '.session' num2str(j) '.ECOGdsrate= fsECOG;']);
eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) '.ECOGtime  = tECOG;']);
catch
eval(['datastruct.k1data.block' num2str(block(i)) '.session' num2str(j) ' = [];']);
end
    end
    end
 
end

