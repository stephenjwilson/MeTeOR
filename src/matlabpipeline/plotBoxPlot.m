function [ ] = plotBoxPlot( data,references,title,root,order )
%plotBoxPlot Plots the box plot data. Not used in figures, but gives a quick readout

if nargin<5
    order=1:length(references); 
end

data=data(:,order);
references=references(order);

h=figure('Position',[100,100,800,800]);
% set(0,'defaultAxesFontName', 'Times');
% set(0,'defaultTextFontName', 'Times');
set(gcf,'visible','off');
hold on;
[~,s2]=size(data);
grp=[];
flat=[];
for i=1:s2
    rel=data(:,i);
    relsize=nnz(rel);
    rel=rel(1:relsize);
    grp=[grp repmat(i,1,relsize)]; %#ok<*AGROW>
    flat=[flat;rel];
end
lw=2;
if mod(s2,2)==0
    totalplots=s2/2;
    c=1;
    for i=1:2:s2
        ax=subplot(1,totalplots,c);

        %ax=subplot('Position',[i/totalplots-1/totalplots+.1,0.1, 0.4, 0.8])
        subset=or(grp==i,grp==i+1);
        %box(ax,'Off');
        bh=boxplot(ax,flat(subset),grp(subset)','Widths',0.9);%,'Labels',references(i:i+1));
        set(bh,'linewidth',lw);
        set(ax,'box','off','color','none')
        set(ax,'fontsize',30);
        if ~isempty(strfind(title,'AUC'))
            m=floor(min(data(data>0))*100)/100-0.1;
            if m>0.4
                m=0.5;
            end
            ylim([0.2,1])
        else
            m=floor(max(data(data>0))*100)/100+.1;
            if m>1
                m=1;
            end
            ylim([0,m])
        end
        c=c+1;
    end
else
    bh=boxplot(flat,grp');
    set(bh,'linewidth',lw);
    set(gca,'box','off','color','none')
    set(gca,'fontsize',30);
    if max(flat)<0.2
        ylim([0,max(flat)+.1])
    else
        ylim([0.2,1])
    end
end
%xlabel('Reference')
%ylabel('Performance')
%ylim([0,1])

%set(gca,'YTick',0:10:120)

hold off;
t = findobj('Tag','Median');
set(t,'Color','k');
t = findobj('Tag','Box');
set(t,'Color','k');
%saveas(h,sprintf('%sBoxPlot_%s_noleg.eps',root,title),'epsc'); 
%set(gca,'XTickMode','auto','XTickLabel',references,'XTick',1:length(references));
saveas(h,sprintf('%sBoxPlot_%s.eps',root,title),'epsc'); 

end

