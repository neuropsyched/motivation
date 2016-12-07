function [mObj,Output,sList] = fbs_trialbypsd(mObj,Output)
% Trial by PSD
%
%
% Stathis Kondylis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tc(1)=tic;
pinfo = mObj.pinfo;
%% Select Recordings
recordstr = cell(size(pinfo));
for i=1:size(pinfo,1)
    recordstr{i,1} = pinfo(i,:).pID;
end
[sList,~] = listdlg('ListString',recordstr,'PromptString','Select recordings to analyze','SelectionMode','multiple');

%% Initialize SPECT
% if nargin<2
%     PSDinfo=struct;
%     save('tempPSD','PSDinfo','-v7.3')
%     n=matfile('tempPSD','writable',true);
%     k=1;
% else
%     k=2;
% end
PSDinfo = struct;
%% Set Conditions

%% RTS
cntrnm = {'Response';'Rsp'};
buffi = 1001:5000;
epoch = {[1601:2000];[2001:2400]};%200 ms before to 400ms afterf0 = [1:1:150]';
%f0 = [1:1:150]';
f0 = 12:23:35;


count=1;
for i = 1:length(sList) %s1
    if pinfo(i,1).Exclude ~= 1
        pID     = pinfo(i,1).pID;
        centers = pinfo(i,1).centers;
        fs      = pinfo(i,1).fs;
        ncntr   = size(centers,2);
        nepoch     = size(epoch,1);
        
        Rsp_str = [pID, '_', 'Rsp'];
        datacell = mObj.(Rsp_str);
        [~,ntrial,nchan] = size(datacell.(Rsp_str)); %ntrial nchannel
        
        PSD = nan(size(f0,2), ntrial, nchan, ncntr, nepoch, 'single'); %[fq, tr, ch, cntr, ep]
        for cntr_i = 1:1 % ncntr
            cntr = centers{cntr_i};
            cntr_str = cell2mat(cntrnm(2,find(strcmp(cntrnm(1,:),cntr))));
            var_str2 = [pID, '_', cntr_str];
            trialdata = datacell.(var_str2);
            
            Nanidx = find(~isnan(trialdata(1,:,1)));
            PSDidx1 = PSD(:,Nanidx,:,cntr_i,1);
            PSDidx2 = PSD(:,Nanidx,:,cntr_i,2);
            
            parfor ch_i = 1:nchan
                chdata = trialdata(buffi,:,ch_i);
                [welchPsd1,~] = pwelch(chdata(epoch{1},Nanidx), 256, 100, f0,fs);
                [welchPsd2,~] = pwelch(chdata(epoch{2},Nanidx), 256, 100, f0,fs);
                PSDidx1(:,:,ch_i) = welchPsd1;
                PSDidx2(:,:,ch_i) = welchPsd2;
            end
            
            PSD(:, Nanidx, :, cntr_i, 1)=PSDidx1;
            PSD(:, Nanidx, :, cntr_i, 2)=PSDidx2;
            
            clc
            fprintf( sprintf('Subject - %s (%d/%d)\n Center - %d/%d\n Time Elapsed - %.2f min\n\n\n', ...
                pID, count, length(sList), cntr_i, ncntr, toc(tc(1))/60  ) )
        end
        
        Output.(pID) = PSD;
        PSDinfo(i,1).pID = pID;
        PSDinfo(i,1).fq  = f0;

        count=count+1;
    end
    % save info
    Output.PSDinfo = PSDinfo;
end
disp('***DONE***')
end


%% SetStimIndex
function index = SetStimInd(cntr,AndConds,OrConds,trial)
for cnd_i=1:size(AndConds,2)
    index{cnd_i}=IndFromCond(AndConds{cnd_i},OrConds{cnd_i},trial);
end
end

function index = FixNanTrials(index,trialdata)
for ind_i=1:size(index,2)
    index{ind_i}(index{ind_i}>size(trialdata,2)) = [];
    index{ind_i}(isnan(trialdata(1,index{ind_i},1)))=[];
end
end

function [P, TStat] = ttest_local (x,y)
%x - [fq, t, tr]
%y - [fq, t]
n      = size(x);
TStat   = zeros(n(1),1,'single');
P       = zeros(n(1),1,'single');

for i=1:n(1)
    x1=squeeze(x(i,:));
    y1=squeeze(y(i,:));
    dim=2;
    
    nx  = size(x1,dim);
    ny  = size(y1,dim);
    s2x = nanvar(x1,[],dim);
    s2y = nanvar(y1,[],dim);
    difference = nanmean(x1,dim) - nanmean(y1,dim);
    dfe = nx + ny - 2;
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    se = sPooled .* sqrt(1./nx + 1./ny);
    ts = difference ./ se;
    p = 2 * tcdf(-abs(ts),dfe);
    %% save into variables
    TStat(i,1)=ts;
    P(i,1)=p;
end

end