function [auc] = GSretro_FastAlt( oldNames,oldNetwork, newNetwork, newNames, predNetwork,predNames)
%GSRetro Does a retrospective analysis on a gold standard.
%   title is the name of the experiment. 
%   All names are cell arrays of strings
%   All networks are matricies
% Assumes symmetrical

%% Figure out mapping to a common set of the first dimension
commonNames=predNames(ismember(predNames,oldNames));
commonNames(strcmp(commonNames,''))=[];

commonNames=newNames(ismember(newNames,commonNames));
commonNames(strcmp(commonNames,''))=[];

[a,old_to_common_mapping]=ismember(commonNames,oldNames);
old_to_common_mapping=old_to_common_mapping(a);

[a,pred_to_common_mapping]=ismember(commonNames,predNames);
pred_to_common_mapping=pred_to_common_mapping(a);

[a,new_to_common_mapping]=ismember(commonNames,newNames);
new_to_common_mapping=new_to_common_mapping(a);

%% Network Manip

predNetwork=predNetwork(pred_to_common_mapping,pred_to_common_mapping);
oldNetwork=oldNetwork(old_to_common_mapping, old_to_common_mapping);
newNetwork=newNetwork(new_to_common_mapping,new_to_common_mapping);

% trim low degree
trimfactor=2;
s1=sum(oldNetwork>0,1);
predNames=predNames(s1>trimfactor);
predNetwork=predNetwork(s1>trimfactor,s1>trimfactor);
oldNetwork=oldNetwork(s1>trimfactor,s1>trimfactor);
newNetwork=newNetwork(s1>trimfactor,s1>trimfactor);
while sum(s1<trimfactor)>0
    s1=sum(oldNetwork>0,1);
    predNames=predNames(s1>trimfactor);
    predNetwork=predNetwork(s1>trimfactor,s1>trimfactor);
    oldNetwork=oldNetwork(s1>trimfactor,s1>trimfactor);
    newNetwork=newNetwork(s1>trimfactor,s1>trimfactor);
    s1=sum(oldNetwork>0,1);
end

nonexistant = tril(oldNetwork==0);
missing_edges = and(tril(newNetwork>0),nonexistant);
nonexistant = (nonexistant - missing_edges)>0;

% Get random indices
n=100000;
num_of_pos = nnz(missing_edges);
num_of_neg = nnz(nonexistant);
pos_rand_indices = randi(num_of_pos, [1,n]);
neg_rand_indices = randi(num_of_neg, [1,n]);

% Isolate Scores
missing_scores = predNetwork(missing_edges);
nonexistant_scores = predNetwork(nonexistant);

% Get vaiables to track auc information
n_prime = sum(missing_scores(pos_rand_indices) > nonexistant_scores(neg_rand_indices));  % Number of times pos score > neg score
n_doubleprime = sum(missing_scores(pos_rand_indices) == nonexistant_scores(neg_rand_indices)); % Number of times pos score = neg score

auc = (n_prime + 0.5 * n_doubleprime) / n ; 
fprintf('Alt AUC: %f\n',full(auc));
end