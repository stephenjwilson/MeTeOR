function [bestks,bestaucs,MCC1s,MCC10s,MCC50s] = runCrossValidation_symmandnonsymm(Networks,AllNetworks,root,ks)

if ~exist(root, 'dir')
  mkdir(root);
end
trimfactor=20;

n=length(ks)*length(Networks);
X=zeros(1000,n);
Y=zeros(1000,n);
AUC=zeros(1,n);

names={};
bestks={0,0,0};
bestaucs={0,0,0};
MCC1s=zeros(length(ks),length(Networks));
MCC10s=zeros(length(ks),length(Networks));
MCC50s=zeros(length(ks),length(Networks));
for o=1:length(Networks)
    for i=1:length(ks)
        %% Get Data
        TestOldNetName=Networks{o}
        data=AllNetworks.(TestOldNetName);
        Network=data{1};
        Names1=data{2};
        Names2=data{2};
        %% Trim Network
        % Trim by number of associations
        s1=sum(Network,1);
        Names1=Names1(s1>trimfactor);
        Network=Network(s1>trimfactor,s1>trimfactor);
        while sum(s1<trimfactor)>0
            s1=sum(Network,1);
            Names1=Names1(s1>trimfactor);
            Network=Network(s1>trimfactor,s1>trimfactor);
            s1=sum(Network,1);
        end
        
        %% Test
        name=sprintf('%s-CrossVal-%s',TestOldNetName,num2str(ks{i}))
        [xs,ys,MCC1,MCC10,MCC50,aucs]=crossvalidation(Network,name,ks{i},1,root);
        MCC1s(i,o)=MCC1;
        MCC10s(i,o)=MCC10;
        MCC50s(i,o)=MCC50;
        X(:,(o-1)*length(ks)+i)=xs;
        Y(:,(o-1)*length(ks)+i)=ys;
        AUC((o-1)*length(ks)+i)=aucs;
        names{(o-1)*length(ks)+i}=name;
        if aucs>bestaucs{o}
            bestaucs{o}=aucs;
            bestks{o}=ks{i};
        end
    end
end
% if length(names)==1
%     names{1}='MeTeOR';
% end

f1=@(name,a) sprintf('%s AUC=%0.2f',name,a);
s1=cellfun(f1, names, num2cell(AUC),'UniformOutput',false);
plot_perfcurve(X,Y,s1,0,'Cross-validation',root);

end