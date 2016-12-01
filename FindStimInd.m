function ind = FindStimInd(AndConds,OrConds,trial,trialdata)

ind = SetStimInd(AndConds,OrConds,trial);
ind = FixNanTrials(ind,trialdata);

end


function ind=SetStimInd(AndConds,OrConds,trial)
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