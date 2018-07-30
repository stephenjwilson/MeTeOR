function [ auc ] = alt_auc( Testnetwork,GS,root,n )
%alt_auc Completes an alternative assessment of auc of the network to a ground truth. 
%   Testnetwork is a network that is being assessed (must be mapped to GS!)
%   GS is a gold standard network (must be mapped to testnetwork!)
if nargin<4
    n=10000; 
end

if ~exist(root, 'dir')
  mkdir(root);
end

positives=and(Testnetwork>0,GS>0);
negatives=and(Testnetwork>0,GS==0);

% Eliminate Duplicates
positives=tril(positives);
negatives=tril(negatives);
% Get Scores
positive_scores = nonzeros(Testnetwork(positives));
negative_scores = nonzeros(Testnetwork(negatives));

% Get random indices
num_of_pos = length(positive_scores);
num_of_neg = length(negative_scores);
pos_rand_indices = randi(num_of_pos, [1,n]);
neg_rand_indices = randi(num_of_neg, [1,n]);

% Get vaiables to track auc information
n_prime = sum(positive_scores(pos_rand_indices) > negative_scores(neg_rand_indices));  % Number of times pos score > neg score
n_doubleprime = sum(positive_scores(pos_rand_indices) == negative_scores(neg_rand_indices)); % Number of times pos score = neg score

auc = (n_prime + 0.5 * n_doubleprime) / n ; 

end

