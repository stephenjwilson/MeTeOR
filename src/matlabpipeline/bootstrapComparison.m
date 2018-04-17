function [ aucs,PRs,PR5s,F1s,F15s ] = bootstrapComparison( Testnetwork,GS,root,title,bootstrapnum )
%bootstrapComparison Completes a bootstrap of the network to a ground truth. The main
%purpose of this is to obtain balanced positive and negative classes so
%that an ROC curve can be accurately calculated
%   Testnetwork is a network that is being assessed (must be mapped to GS!)
%   GS is a gold standard network (must be mapped to testnetwork!)
%   bootstrap num is the num of bootstraps to try and obtain. If there are
%   too few positives, the bootstrap number will be adjusted.
if nargin<5
    bootstrapnum=100; 
end

if ~exist(root, 'dir')
  mkdir(root);
end

aucs=zeros(bootstrapnum,1);
PRs=zeros(bootstrapnum,1);
PR5s=zeros(bootstrapnum,1);
F1s=zeros(bootstrapnum,1);
F15s=zeros(bootstrapnum,1);

positives=and(Testnetwork>0,GS>0);
negatives=and(Testnetwork>0,GS==0);

fprintf('Positives: %i. Negatives: %i\n',nnz(positives),nnz(negatives));

sizeofboot=floor(nnz(positives)/bootstrapnum);
sizeofbootneg=floor(nnz(negatives)/sizeofboot);
if sizeofboot<100
    disp('Iteration number is too high for the number of positives. Adjusting...');
end
while sizeofboot<100
    bootstrapnum=floor(bootstrapnum/2);
    sizeofboot=floor(nnz(positives)/bootstrapnum);
    sizeofbootneg=floor(nnz(negatives)/sizeofboot);
end
aucs=aucs(1:bootstrapnum);
PRs=PRs(1:bootstrapnum);
PR5s=PR5s(1:bootstrapnum);
F1s=F1s(1:bootstrapnum);
F15s=F15s(1:bootstrapnum);

cvClassespos=floor(rand(nnz(positives),1)*bootstrapnum)+1;
cvClassesneg=floor(rand(nnz(negatives),1)*sizeofbootneg)+1;

labelspos=Testnetwork(positives);
labelsneg=Testnetwork(negatives);
for i=1:bootstrapnum
    rankpos=labelspos(cvClassespos==i);
    rankneg=labelsneg(cvClassesneg==i);
    ranking=[rankpos;rankneg];
    classes=[ones(length(rankpos),1); zeros(length(rankneg),1)];
    try
        [~,~,~,auc_ROC]=perfcurve(classes,ranking,1);    
        [x_precRecall,y_precRecall,~,~]=perfcurve(classes,ranking,1,'XCrit','TPR','YCrit','PPV');
    catch
        auc_ROC=0.5;
        x_precRecall=[0];
        y_precRecall=[1];
    end
    aucs(i)=auc_ROC;
    
    r=0.1;
    ind=find(x_precRecall>r,1);
    if isempty(ind)
        p=0;
        r=1;
    else
        p=y_precRecall(ind);
        r=x_precRecall(ind);
    end
    f1=2*(p*r)/(p+r);
    PRs(i)=p;
    F1s(i)=f1;
    r=0.5;
    ind=find(x_precRecall>r,1);
    if isempty(ind)
        p=0;
        r=1;
    else
        p=y_precRecall(ind);
        r=x_precRecall(ind);
    end
    f1=2*(p*r)/(p+r);
    PR5s(i)=p;
    F15s(i)=f1;
end
h=figure('Position',[100,100,800,800]);
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(gcf,'visible','off');
set(gca,'fontsize',25);
hold on;
boxplot([aucs,PRs,F1s,PR5s,F15s],{'AUC','PR(0.1)','F1(0.1)','PR(0.5)','F1(0.5)'})
xlabel('Metric')
ylabel('Performance')
ylim([0,1])
hold off;
saveas(h,sprintf('%sBoxPlot_%s.eps',root,title),'epsc'); 
end

