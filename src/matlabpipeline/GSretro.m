function [BAUCs,BPR,BPR5,BF1,BF15 ] = GSretro( oldNames1,oldNames2,predNetwork,predNames1,predNames2, GSoldNetwork, GSoldNames1, GSoldNames2, GSnewNetwork, GSnewNames1,GSnewNames2, title,root)
%GSRetro Does a retrospective analysis on a gold standard.
%   title is the name of the experiment. 
%   All names are cell arrays of strings
%   All networks are matricies

%% Figure out mapping to a common set of the first dimension
commonNames1=predNames1(ismember(predNames1,GSoldNames1));
commonNames1(strcmp(commonNames1,''))=[];

commonNames1=GSnewNames1(ismember(GSnewNames1,commonNames1));
commonNames1(strcmp(commonNames1,''))=[];

commonNames1=oldNames1(ismember(oldNames1,commonNames1));
commonNames1(strcmp(commonNames1,''))=[];

[a,old_to_common_mapping1]=ismember(commonNames1,oldNames1);
old_to_common_mapping1=old_to_common_mapping1(a);

[a,pred_to_common_mapping1]=ismember(commonNames1,predNames1);
pred_to_common_mapping1=pred_to_common_mapping1(a);

[a,GSold_to_common_mapping1]=ismember(commonNames1,GSoldNames1);
GSold_to_common_mapping1=GSold_to_common_mapping1(a);

[a,GSnew_to_common_mapping1]=ismember(commonNames1,GSnewNames1);
GSnew_to_common_mapping1=GSnew_to_common_mapping1(a);

%% All of the second dimension
commonNames2=predNames2(ismember(predNames2,GSoldNames2));
commonNames2(strcmp(commonNames2,''))=[];

commonNames2=GSnewNames2(ismember(GSnewNames2,commonNames2));
commonNames2(strcmp(commonNames2,''))=[];

commonNames2=oldNames2(ismember(oldNames2,commonNames2));
commonNames2(strcmp(commonNames2,''))=[];

[a,pred_to_common_mapping2]=ismember(commonNames2,predNames2);
pred_to_common_mapping2=pred_to_common_mapping2(a);

[a,GSold_to_common_mapping2]=ismember(commonNames2,GSoldNames2);
GSold_to_common_mapping2=GSold_to_common_mapping2(a);

[a,GSnew_to_common_mapping2]=ismember(commonNames2,GSnewNames2);
GSnew_to_common_mapping2=GSnew_to_common_mapping2(a);

disp('Found Mappings');

%% Network Manip
% Get GS
GSnewNetwork=GSnewNetwork(GSnew_to_common_mapping1,GSnew_to_common_mapping2)>0;
GSoldNetwork=GSoldNetwork(GSold_to_common_mapping1,GSold_to_common_mapping2)>0;
% Get pred
predNetwork=predNetwork(pred_to_common_mapping1,pred_to_common_mapping2);
% Use GSold to eliminate things from GSnew and from pred
mask=GSoldNetwork>0;
GSnewNetwork(mask)=-1;
predNetwork(mask)=-1;
clear mask;
positives=and(predNetwork>0,GSnewNetwork>0);
negatives=and(predNetwork>0,GSnewNetwork==0);
auc_alt = alt_auc(predNetwork,GSnewNetwork,root,10000);
fprintf('Alt AUC: %f\n',auc_alt);
%% Bootstrap
[aucs,PRs,PR5s,F1s,F15s ]=bootstrapComparison(predNetwork,GSnewNetwork,strcat(root,'/BoxPlot/'),title);

BAUCs=aucs;
BPR=PRs;
BPR5=PR5s;
BF1=F1s;
BF15=F15s;

%% Continue with other ROC / PR
fprintf('Positives: %f\n',nnz(positives));
fprintf('Negatives: %f\n',nnz(negatives));

if length(GSoldNames1)==length(GSoldNames2)
    if sum(~strcmp(GSoldNames1,GSoldNames2))==0 %%%%%
        positives=tril(positives);
        negatives=tril(negatives);
    end
end
rankpos = predNetwork(positives);
rankneg = predNetwork(negatives);
ranking=[rankpos;rankneg];
classes=[ones(length(rankpos),1); zeros(length(rankneg),1)];
%% Get ROC curve
[x_ROC,y_ROC,~,auc_ROC]=perfcurve(classes,ranking,1);
disp(auc_ROC)

%% Reduce dimensionality if needed
if length(x_ROC)>1000
    x_ROC=reduceVector(x_ROC,1000);
    y_ROC=reduceVector(y_ROC,1000);
else
    x_ROC(end:1000)=1;
    y_ROC(end:1000)=1;
end


%% Precision recall curves
[x_precRecall,y_precRecall,~,~]=perfcurve(classes,ranking, 1, 'XCrit','TPR','YCrit','PPV');

if length(x_precRecall)>1000
    x_precRecall=reduceVector(x_precRecall,1000);
    y_precRecall=reduceVector(y_precRecall,1000);
else
    x_precRecall(end:1000)=x_precRecall(end);
    y_precRecall(end:1000)=y_precRecall(end);
end

%plot curves
plot_perfcurve(x_ROC,y_ROC,{sprintf('%s AUC=%0.2f',title, auc_ROC)},0,title,root);
plot_perfcurve(x_precRecall,y_precRecall,{title},1,title,root);


%% Output top predictions
pr=predNetwork(positives);
[row,col,v]=find(positives);
l1=commonNames1(row);
l2=commonNames2(col);
data=table(l1,l2,pr);
data.Properties.VariableNames = {'Entity1','Entity2','Weight'};
data = sortrows(data,'Weight','descend');
writetable(data,sprintf('%sPredlistRetro_%s.txt',root,title),'Delimiter','\t');

s2={sprintf('%s',title)};
%% MetaBoxplot
boxroot=strcat(root,'/BoxPlot/');
plotBoxPlot(BAUCs,s2,strcat(title,'AUC'),boxroot);
plotBoxPlot(BPR,s2,strcat(title,'PR1'),boxroot);
plotBoxPlot(BPR5,s2,strcat(title,'PR5'),boxroot);
plotBoxPlot(BF1,s2,strcat(title,'F11'),boxroot);
plotBoxPlot(BF15,s2,strcat(title,'F15'),boxroot);

save(sprintf('GSRetro_%s.mat', title),'BAUCs','BPR','BPR5','BF1','BF15','x_ROC','y_ROC','x_precRecall', 'y_precRecall','auc_alt')

end

