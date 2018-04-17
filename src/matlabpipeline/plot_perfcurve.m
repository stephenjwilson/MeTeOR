function [ ] = plot_perfcurve(X,Y,labels,prec,title,root,order,suborder,colors)
%plot_ROC Plots a performance curve
%   ROC curve if prec=0, precision recall if prec=1. Works on up to 4 lines
%   of a single type
    %% Set Up
    if nargin<6
        root='./';
    elseif nargin<7
        order=[];
    elseif nargin<8
        suborder=[];
    end
    if nargin<9
        colors={[1,0,0],[0,0,1],[0.2,0.5,0.2],[0.5,0.5,0.5]};
    end
    h=figure('Position',[100,100,800,800]);
    set(0,'defaultAxesFontName', 'Times');
    set(0,'defaultTextFontName', 'Times');
    set(gcf,'visible','off');
    %set(gcf,'defaultAxesLineStyleOrder','-|:');
    set(gca,'fontsize',20);
    

    
    %% Plot Data
    hold on;
    if isempty(order)
        plot(X,Y,'LineWidth',2)%,'.');% X and Y must be matricies that are equal sizes and
             % columns will be plotted
    else
        ind=zeros(length(labels),1);
        u=unique(order);
%         colors={'r','b','g','c','m','y'};
        
        dashes={'-','--',':','-.'};
        count=1;
        for i=1:length(u)
            x=X(:,order==u(i));
            y=Y(:,order==u(i));
            ls=labels(order==u(i));
            if ~isempty(suborder)
                sorder=suborder(order==u(i));
                
                if sorder(1)==2
                    
                    x(:,[1,2])=x(:,[2,1]);
                    y(:,[1,2])=y(:,[2,1]);
                    ls([1,2])=ls([2,1]);
                end
            end
            [~,w]=size(x);
            
            for o=1:w
                 plot(x(:,o),y(:,o),'LineWidth',2,'Color',colors{i},'LineStyle',dashes{o});%,'.');% X and Y must be matricies that are equal sizes and
            end
            
            %Really messy label fixing
            for o=1:length(ls)
                tmp=find(strcmp(labels,ls{o}));
                if length(tmp)>1
                    ind(count)=tmp(i);
                else
                    ind(count)=tmp;
                end
                count=count+1;
            end
        end
        labels=strrep(labels(ind),'_','');
    end

    %title(name,'FontSize',16);
    
    if prec==1
        xlabel('Recall','FontSize',30);
        axis([0 1 0 1]);
        ylabel('Precision','FontSize',30);
        saveas(h,sprintf('%sPR_%s_noleg.eps',root,title),'epsc');
        leg=legend(labels,'location','North');
        set(leg,'FontSize',20);
        saveas(h,sprintf('%sPR_%s.eps',root,title),'epsc');
        % Save Data
        dlmwrite(sprintf('%sPR%s_data.out',root,title),[X Y], '\t');
        %saveas(h,sprintf('%sPR_%s.fig',root,title));
    elseif prec==2
        xlabel('Recall','FontSize',30);
        %axis([0 1 0 1]);
        ylabel('Accuracy','FontSize',30);
%         leg=legend(labels,'location','North');
%         set(leg,'FontSize',20);
        saveas(h,sprintf('%sacc_%s.eps',root,title),'epsc');
        
        %saveas(h,sprintf('%sacc_%s.fig',root,title));
    else
        plot([0,1],[0,1],'k');
        labels{length(labels)+1}='Random';
        
        xlabel('False Positive Rate','FontSize',30);
        ylabel('True Positive Rate','FontSize',30);    
        saveas(h,sprintf('%sROC_%s_noleg.eps',root,title),'epsc');
        leg=legend(labels,'location','SouthEast');
        set(leg,'FontSize',20);        
        % Save Data
        dlmwrite(sprintf('%sROC%s_data.out',root,title),[X Y], '\t');
        saveas(h,sprintf('%sROC_%s.eps',root,title),'epsc'); 
%         saveas(h,sprintf('%sROC_%s.fig',root,title));
    end
   %set(groot,'defaultAxesLineStyleOrder','remove')
    %set(groot,'defaultAxesColorOrder','remove')
end

