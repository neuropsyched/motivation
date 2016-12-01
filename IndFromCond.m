function [ind] = IndFromCond (AndCond, OrCond, trial)

ind=[];
if isempty(OrCond)
    ind=union(trial.Correct,trial.Error);
else
    for i=1:size(OrCond,2)
        ind=union(ind,trial.(OrCond{i}));
    end
end

for i=1:size(AndCond,2)
    ind=intersect(ind,trial.(AndCond{i}));
end

[r,c]=size(ind);
if c>r, ind=ind'; end