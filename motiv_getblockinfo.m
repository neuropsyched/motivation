function [info] = motiv_getblockinfo(info)

% datfile refers to the file name MT0001 etc...
% blockname refers to the block name 

if ~isfield(info,'block')
datfiles = zeros(1,length(info.reply)); 
for i=1:length(info.reply)
    datfiles(i) = str2double(info.reply(i));
end
clear i
end

nev = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.nev')));
ns2 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns2')));
ns5 = info.TrellisFilenames(~cellfun(@isempty,strfind(info.TrellisFilenames,'.ns5')));

% Find Start and Stop Times
if ~isfield(info,'tstart') || ~isfield(info,'tstop')
for i = 1:length(datfiles)
    if ~(length(nev)<i)
        [EventTimes] = GetEventData(fullfile(info.trellispath,nev{i}));
        % Plot channels
        if ~isempty(EventTimes) % || ~EventTimes==0
            disp(['Displaying EventTimes for datafile ' info.reply{i}])
            h = figure; plot(EventTimes*[1 1], ylim, 'k')
            title(['datfile' num2str(datfiles(i))])
            waitfor(h) % pause while figure is open
            disp('Choose start and stop times')
            blocknames = strsplit(cell2mat(inputdlg('Enter block names')),','); % :)
            tstart = strsplit(cell2mat(inputdlg('Enter start time')),','); % :)
            tstop = strsplit(cell2mat(inputdlg('Enter stop time')),',');
            % str2num(char(cellfun(@(x) (x(1:strfind(tstart{1},',')-1)), tstart,'Uniformoutput',false)))
        else
            disp(['Empty EventTimes for datafile ' info.reply{i}])
            blocknames = 'na';
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
    
    info.datfiles = datfiles;
    info.blocknames{i} = blocknames;
    info.tstart{i} = tstart;
    info.tstop{i} = tstop;   
end
save(fullfile(info.patientdir,'motiv_info.mat'),'info')
else
end

