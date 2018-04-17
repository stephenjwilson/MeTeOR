function [ xs,ys,MCC1,MCC10,MCC50,aucs ] = crossvalidation( testnetwork,title,k,symm,root,crossValidationFold)
%CROSSVALIDATION Modified from Ben Bachman's Crossvalidation code
if nargin<6
    crossValidationFold=10;
end
if nargin<4
    symm=1;
end
if nargin<5
    root='./';
end
%ks=[1 5 10 20 50 100 200 500 1000 2000];

%assign each edge to a crossvalidation class, make it symmetric
cvClasses=floor(rand(size(testnetwork))*crossValidationFold)+1;
% %ignore nodes with fewer than 20 edges, they'll be used to build the MF
% % model, but not tested
removeThese=sum(testnetwork>0)<20;
cvClasses(removeThese,:)=-1;
cvClasses(:,removeThese)=-1;
% % make it so we remove both directions at the same time
% cvClasses=triu(cvClasses)+triu(cvClasses)';


%remove diagonals, none of these in string, 
if symm
    cvClasses(eye(size(cvClasses))==1)=-1;
end

scores=zeros(size(testnetwork));
classes=testnetwork>0;

xs=zeros(1000,1);
ys=zeros(1000,1);

aucs=zeros(1,1);
sum(sum(testnetwork<0));
% disp('Trying nnmf')
% nnmf(testnetwork,200)
try
    load([root num2str(crossValidationFold) '_crossvalidation_results_k' num2str(k) '_' title '.mat']);
catch
    for cvClass=[1:crossValidationFold]
        cvClass
        network_copy=testnetwork;
        testset=(cvClasses==cvClass);
        %hide testset
        network_copy(testset)=0;
%         network_copy<0;
%         sum(sum(network_copy<0));
        seed = 1;
        n = 10;
        opt = statset('MaxIter',10,'Display','final','Streams',RandStream.create('mrg32k3a','NumStreams',n,'Seed',seed),'UseSubstreams',true);
        [W,H] = nnmf(network_copy,k,'algorithm','mult','options',opt);
%         [W,H]=nnmf(network_copy,k);
        result=W*H;
        %average it out
        if symm
            result=(result+result')/2;
        end
        %store results
        scores(testset)=result(testset);
    end
    %make it all a big vector for perfcurve
    if symm
        scores=reshape(scores,1,size(testnetwork,1)^2);
    else
        [w,h]=size(testnetwork);
        scores=reshape(scores,1,w*h);
    end

    clear network_copy W H result testset
    %% Calculate success
    % Calculate AUC
    [xs,ys,t,aucs]=perfcurve(classes(cvClasses~=-1),scores(cvClasses~=-1),1);
    ind1=find(xs>0.01,1); % 10% Recall
    ind10=find(xs>.1,1); % 10% Recall
    ind50=find(xs>.5,1);
    xs=reduceVector(xs,1000);
    ys=reduceVector(ys,1000);
    disp('Done with AUC') 
    
    % Calculate MCC at vals
    labels=classes(cvClasses~=-1)';
    pred1=scores(cvClasses~=-1)>t(ind1);
    pred10=scores(cvClasses~=-1)>t(ind10);
    pred50=scores(cvClasses~=-1)>t(ind50);
    MCC1=calcMCC(labels,pred1);
    MCC10=calcMCC(labels,pred10);
    MCC50=calcMCC(labels,pred50);

    save([root num2str(crossValidationFold) '_crossvalidation_results_k' num2str(k) '_' title '.mat']);
end
clear cvClasses
clear scores
    
end
