function [matObj,testPSDout,sList] = fbs_testpsd(matObj,PSD,sList)
% Plot PSD
%
%
% Stathis Kondylis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pinfo    =  matObj.pinfo;
PSDinfo  =  PSD.PSDinfo;
%% Select Recordings in fbs_trialbypsd and pass sList %otherwise choose
if nargin==2
    recordstr = cell(size(pinfo));
    for i=1:size(pinfo,1)
        recordstr{i,1} = pinfo(i,:).pID;
    end
    [sList,~] = listdlg('ListString',recordstr,'PromptString','Select recordings to analyze','SelectionMode','multiple');
elseif nargin==1
    error('not enough inputs')
end

%%
AndConds{1}={'Correct','Right'};
    OrConds{1}={'Win','Squeeze','Go','GoFast'};
AndConds{2}={'Correct','Left'};
    OrConds{2}={'Win','Squeeze','Go','GoFast'};

col = [1,0,0;0,0,1;0,1,0];

Bands = [8,32;76,100];
nbnd = size(Bands,1);
bndN = {'LFB','HFB'}

clc
for i = 1:length(PSDinfo)
    if pinfo(i,1).Exclude ~= 1
        
        pID   = PSDinfo(i,1).pID;
        fq    = PSDinfo(i,1).fq;
        trial = pinfo(i,1).trial;
        psdnow     = PSD.(pID);            %[fq, tr, ch, cntr, ep]
        labels = pinfo(i,1).labels;
        %labels = pinfo(sbj_i,1).labelsBP.labels2;
        
        Rsp_str = [pID, '_', 'Rsp'];
        datacell = matObj.(Rsp_str);
        trialdata = size(datacell.(Rsp_str)); %ntrial nchannel
        
        cntr='Response';
        ind=SetStimInd(cntr,AndConds,OrConds,trial);%Only Need This When Calculating ZMI
        ind=FixNanTrials(ind,trialdata);
        ncnd = size(ind,2);
        
        [~, ntr, nch, ncntr, nep]=size(psdnow);
        
        nr = ceil(sqrt(nch));
        nc = nr;
        c1i=1:nc:nr*nc;
        rlasti=nch-nr+1:1:nch;
        
        Stat = zeros(nch,ncnd,nbnd,2,'single'); %[channels, conditions, band, P/TS, epoch]
        
        for bnd_i = 1:nbnd
            fh(i)=figure(1); 
            %set(gcf,'windowstyle','docked'); 
            clf
            sh{i}=gobjects(nch,1);
            
            for ch_i = 1:nch
                
                sh{i}(ch_i,1) = subplot(nr,nc,ch_i); hold on
                %sh{sbj_i}(ch_i,2) = axes; hold on
                %pos=get(sh{sbj_i}(ch_i,1),'position')
                
                for cnd_i = 1:ncnd
                    
                    Arts1 = pinfo(i,1).Arts{ch_i,1};
                    Arts2 = pinfo(i,1).Arts{ch_i,2};
                    Arts = union(Arts1,Arts2);
                    
                    ind1 = setdiff(ind{cnd_i},Arts);
                    ind2 = setdiff(ind{cnd_i},Arts);
                    
                    if ~isempty(ind1) && ~isempty(ind2)
                        P1 = squeeze(psdnow(:,ind1,ch_i,1,1));
                        P2 = squeeze(psdnow(:,ind2,ch_i,1,2));
                        P3 = squeeze(psdnow(:,ind2,ch_i,1,3));
                        P11 = P1( Bands(bnd_i,1):Bands(bnd_i,2) ,: );
                        P22 = P2( Bands(bnd_i,1):Bands(bnd_i,2) ,: );
                        P33 = P3( Bands(bnd_i,1):Bands(bnd_i,2) ,: );
                        [A1] = KMtest_local ( P11, P33 );
                        [p1, ts1] = ttest_local (A1,bnd_i);
                        [A2] = KMtest_local ( P22, P33 );
                        [p2, ts2] = ttest_local (A2,bnd_i);
                        
                        Stat(ch_i,cnd_i,bnd_i,1,1) = p1;
                        Stat(ch_i,cnd_i,bnd_i,2,1) = ts1;
                        Stat(ch_i,cnd_i,bnd_i,1,2) = p2;
                        Stat(ch_i,cnd_i,bnd_i,2,2) = ts2; %[ch,cnd,bnd,p/ts,ep]
                        
                        %if p < (0.05), mrk = '*'; else mrk = '.'; end%Bonferroni corrected threshold
                        %ph1 = plot(sh{sbj_i}(ch_i,1), cnd_i, ts, 'color', col(cnd_i,:), 'marker',mrk,'markersize',15);
                    end
                end
%                 title(labels{ch_i})
%                 if ~any(ismember(c1i,ch_i)), set(gca,'ytick',[]), end
%                 if ~any(ismember(rlasti,ch_i)), set(gca,'xtick',[]), end
            end
            
%             set(sh{sbj_i},'xlim',[0,3],'ylim',[-10,10])
%             suptitle(sprintf('Subject: %s\nConditions: Right, Left\nBand: %s',pID, bndN{bnd_i}))
%             
%             
            
            svn = [pID,'_stat'];
            testPSDout.(svn) = Stat;
            
%             ft=sprintf('%s_%s_PSD',pinfo(sbj_i).pID,bndN{bnd_i});
%             ft1=fullfile(cd,'PSD Figures',ft);
%             %FOR SAVE
%             export_fig(ft1,'-bmp','-r50')
            
        end
    end
    disp(sprintf('Subject %d',i))
end

end




%%
    
function ind=SetStimInd(cntr,AndConds,OrConds,trial)
for cnd_i=1:size(AndConds,2)
    ind{cnd_i}=IndFromCond(AndConds{cnd_i},OrConds{cnd_i},trial);
end
end

function ind=FixNanTrials(ind,trialdata)
for ind_i=1:size(ind,2)
    ind{ind_i}(ind{ind_i}>size(trialdata,2)) = [];
    ind{ind_i}(isnan(trialdata(1,ind{ind_i},1)))=[];
end
end

function [p, tval] = ttest_local (x,tail)

m=0;
dim=2;
nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
df = max(samplesize - 1,0);
xmean = nanmean(x,dim);
sdpop = nanstd(x,[],dim);
ser = sdpop ./ sqrt(samplesize);
tval = (xmean - m) ./ ser;

if tail == 1
    p = tcdf(tval, df);%Lt 1 tail
elseif tail == 2
    
    p = tcdf(-tval, df);%Rt 1 tail
elseif tail==3
    p = 2 * tcdf(-abs(tval), df);%2 tail
end



end


function [A] = KMtest_local (x,y)

mxtr = max(size(x,2),size(y,2));

mnx = nanmean(x,1);
mny = nanmean(y,1);
xystd = nanstd([x;y],[],1);

Nx = size(x,1);
Ny = size(y,1);

A = ( ( (mnx - mny).^3 ) ./ ( abs(mnx - mny) .* xystd.^2  ) ) * ( Nx*Ny/(Nx+Ny)^2  );

end
