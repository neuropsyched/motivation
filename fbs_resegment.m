sbj_list = {'RT','VK','JD','MM','WO','SE','MK','JA','JM'};
flt_fq = [13,20];
varmod = 'B1';


for sbj_i = [4,5]


load(fullfile( 'MV2_data', sprintf('%s_noref', sbj_list{sbj_i} ) ) );
load(fullfile( 'MV2_Trialdata', sprintf('%s_trial', sbj_list{sbj_i} ) ) );



%% Segment

s=1:size(referenced,2);
if ~isempty(varmod)
LFP = [eegfilt([referenced(:,s)]',1000,flt_fq(1),flt_fq(2))]';
elseif isempty(varmod)
LFP = referenced(:,s);
end

clear referenced
labels=labels(s);

trialdata={};
trialforce={};


%-------------------% Set Stimulus Times------------------------
centers = {'Cue','Command','Response','Cue'};
tcntr = {[-1.2*fs+1:2.2*fs], [-1.2*fs+1:2.2*fs], [-1.2*fs+1:1.2*fs], [-.55*fs+1:0.2*fs]};
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
    trialdata{cntr_i}=nan(length(tt),length(stim),length(labels));
    trialforce{cntr_i}=nan(length(tt),length(stim),2);
    
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

svname=fullfile('MV2_Trialdata', [sbj_list{sbj_i} '_trial']);

if ~isempty(varmod)
eval(sprintf('%s_trialdata = trialdata;',varmod));
eval(sprintf('%s_tcntr = tcntr;',varmod));
save(svname,'-append',sprintf('%s_trialdata',varmod))
elseif isempty(varmod)
    save(svname,'-append','trialdata','tcntr','trialforce')
end



clc
fprintf(sprintf('Completed - %s\n     %d / %d', sbj_list{sbj_i}, sbj_i, 9));

end