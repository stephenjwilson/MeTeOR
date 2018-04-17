function [MCC] = calcMCC(labels,pred)
%MCC Calculates the Matthews Correlation Coefficient
%   Detailed explanation goes here
    [~,cm,~,~] = confusion(labels,pred);
    
    TP=cm(2,2); %sum(pred(labels>0))
    FP=cm(1,2); %sum(pred(labels==0))
    FN=cm(2,1); %sum(labels(pred==0)==1)
    TN=cm(1,1); %sum(labels(pred==0)==0)
    MCC=((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end

