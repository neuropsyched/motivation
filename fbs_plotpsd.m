function fbs_plotpsd(matObj,PSD,sList)
% Plot PSD
%
%
% Stathis Kondylis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pinfo    =  matObj.pinfo;
PSDinfo  =  PSD.PSDinfo;
%% Select Recordings
if ~exist('sList')
recordstr = cell(size(pinfo));
for i=1:size(pinfo,1)
    recordstr{i,1} = pinfo(i,:).pID;
end
[sList,~] = listdlg('ListString',recordstr,'PromptString','Select recordings to analyze','SelectionMode','multiple');
end
%%
AndConds{1}={'Right','Correct'};
OrConds{1}={'Win','Squeeze','Go','GoFast'};
AndConds{2}={'Left','Correct'};
OrConds{2}={'Win','Squeeze','Go','GoFast'};
col = [1,0,0;0,0,1];

for i = 1:length(sList)
    sbj = sList(i); % Subject index
    figure; set(gcf,'windowstyle','normal'); clf
    pID   = pinfo(sbj,1).pID;
    fq    = PSDinfo(sbj,1).fq;
    trial = pinfo(sbj,1).trial;
    psdnow     = PSD.(pID);           %[fq, tr, ch, cntr, ep]
    labels = pinfo(sbj,1).labels;
    
    Rsp_str = [pID, '_', 'Rsp'];
    trialdata = matObj.(Rsp_str);
    if isstruct(trialdata)            %catch if data is double nested
        trialdata = trialdata.(Rsp_str);
    end
    
    cntr = 'Response';
    %    ind=SetStimInd(cntr,AndConds,OrConds,trial);%Only Need This When Calculating ZMI
    %    ind=FixNanTrials(ind,trialdata);
    index = FindStimInd(AndConds,OrConds,trial,trialdata);
    ncnd = size(index,2);
    
    [~, ~, nch, ~, ~] = size(psdnow);
    %  [fq, tr, ch, cntr, ep]
    nr = ceil(sqrt(nch));
    nc = nr;
    c1i=1:nc:nr*nc;
    rlasti=nch-nr+1:1:nch;
    
    sh{i}=gobjects(nch,1);
    
    for ch_i = 1:nch
        
        sh{i}(ch_i,1) = subplot(nr,nc,ch_i); hold on
        %sh{sbj_i}(ch_i,2) = axes; hold on
        %pos=get(sh{sbj_i}(ch_i,1),'position')
        
        for cnd_i = 1:ncnd
            
            Arts1 = pinfo(i,1).Arts{ch_i,1};
            Arts2 = pinfo(i,1).Arts{ch_i,2};
            
            ind1 = setdiff(index{cnd_i},Arts1);
            ind2 = setdiff(index{cnd_i},Arts2);
            
            P1 = squeeze(psdnow(:,ind1,ch_i,1,1));
            P2 = squeeze(psdnow(:,ind2,ch_i,1,2));
            %[p,ts] = ttest_local ( P2, P1 );
            [A] = KMtest_local ( P2, P1 );
            
            ph1 = plot(sh{i}(ch_i,1), fq, A, 'color', col(cnd_i,:));
            %ph2 = plot(sh{sbj_i}(ch_i,2), fq, p, 'color', col(cnd_i,:), 'linestyle', '--');
            %set(sh{sbj_i}(ch_i,:),'position',pos)
            %set(sh{sbj_i}(ch_i,2),'color','none','xtick',[],'yaxislocation','right')
            %axes(sh{sbj_i}(ch_i,2))
            
        end
        title(labels{ch_i})
        if ~any(ismember(c1i,ch_i)), set(gca,'ytick',[]), end
        if ~any(ismember(rlasti,ch_i)), set(gca,'xtick',[]), end
    end
    
    set(sh{i},'xlim',[0,max(fq)],'ylim',[-.5,.5])
    suptitle(pID)
    
    ft=sprintf('%s_PSD',pinfo(i).pID);
    ft1=fullfile(cd,'PSD Figures',ft);
    %FOR SAVE
%     export_fig(ft1,'-bmp','-r50')
    
end
end


%%
function index = SetStimInd(cntr,AndConds,OrConds,trial)
for cnd_i=1:size(AndConds,2)
    index{cnd_i} = IndFromCond(AndConds{cnd_i},OrConds{cnd_i},trial);
end
end

function ind = FixNanTrials(ind,trialdata)
for ind_i=1:size(ind,2)
    ind{ind_i}(isnan(trialdata(1,ind{ind_i},1)))=[];
end
end

function [P, TS] = ttest_local (x,y)
%x - [fq, t, tr]
%y - [fq, t]
sz=size(x);
TS=zeros(sz(1),1,'single');
P=zeros(sz(1),1,'single');

for i=1:sz(1)
    x1=squeeze(x(i,:));
    y1=squeeze(y(i,:));
    dim=2;
    
    nx=size(x1,dim);
    ny=size(y1,dim);
    s2x = nanvar(x1,[],dim);
    s2y = nanvar(y1,[],dim);
    difference = nanmean(x1,dim) - nanmean(y1,dim);
    dfe = nx + ny - 2;
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    se = sPooled .* sqrt(1./nx + 1./ny);
    ts = difference ./ se;
    p = 2 * tcdf(-abs(ts),dfe);
    
    TS(i,1)=ts;
    P(i,1)=p;
end

end


function [A] = KMtest_local (x,y)

mnx = nanmean(x,2);
mny = nanmean(y,2);
xystd = nanstd([x,y],[],2);
Nx = size(x,2);
Ny = size(y,2);

A = ( ( (mnx - mny).^3 ) ./ ( abs(mnx - mny) .* xystd.^2  ) ) * ( Nx*Ny/(Nx+Ny)^2  );



end
