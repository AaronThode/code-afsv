function [N,printname,Ibin,hout,hprint]=hist2D(X,axis1,axis2,feature_names,contour_chc,scale_chc,trim_chc,no_margin)
%function [N,printname,Ibin]=hist2D(X,axis1,axis2,feature_names,contour_chc,scale_chc,trim_chc)
%
%  Input:
%%% X: [2 Npoint] data sample;
%%%  axis1, axis2: desired histogram bins..
%%%  feature_names: cell array of names of featurs for plotting
%%%  contour_chc: if 1, image, if 2, contourf
%%%  scale_chc:  if 1, linear, if 2, log10.  Default is log10
%%%  trim_chc:   if 1 don't remove data outside axis1 or axis2.  If 2 remove outlier data from  histograms.
%%%  no_margin:  if true don't plot marginal distributions.
%%%
%%% Output:
%%% N : matrix of binned counts; rows with axis1, columns with axis2
%%%      Note that 'x' axis is axis2 feature, etc.
%%% Ibin:  [2 Npoint] matrix that gives [I; J] indicie where sample was placed in N
%%%          First row associated with rows of N (axis 1); second row associated with columns of N
%%% hout: handles to three axes...
%%% hprint: handle to figure
if size(X,1)>2
    error('X needs to have only two rows');
end
Ibin=zeros(size(X));

if ~exist('contour_chc')|isempty(contour_chc)
    contour_chc=1;
end
if ~exist('scale_chc')|isempty(scale_chc)
    scale_chc=2;
end

if ~exist('trim_chc')|isempty(trim_chc)
    trim_chc=1;
end

if ~exist('no_margin')|isempty(no_margin)
    no_margin=false;
end

if trim_chc==2
    Igood=find(X(1,:)>=axis1(1)&X(1,:)<=axis1(end)&X(2,:)>=axis2(1)&X(2,:)<=axis2(end));
    X=X(:,Igood);
    
end

for I=1:length(axis1)
    
    if I==1
        Igood=find(X(1,:)<=axis1(1));
    elseif I==length(axis1);
        Igood=find(X(1,:)>axis1(end));
    else
        Igood=find(X(1,:)>axis1(I)&X(1,:)<=axis1(I+1));
    end
    
    for J=1:length(axis2)
        
        if J==1
            Igood2=find(X(2,Igood)<=axis2(1));
        elseif J==length(axis2);
            Igood2=find(X(2,Igood)>axis2(end));
        else
            Igood2=find(X(2,Igood)>axis2(J)&X(2,Igood)<=axis2(J+1));
        end
        N(I,J)=length(Igood2);
        Icell=Igood(Igood2);  %indicies of X that fill cell I,J
        if ~isempty(Icell)
            Ibin(:,Icell)=[I; J]*ones(1,length(Icell));
        end
    end
    
end


hprint=figure;

%%%%%%%%%%%%%%%%Margin 1
if ~no_margin
    subplot(2,2,1);
    [N1,X1]=hist(X(1,:),axis1);
    barh((X1),(N1)/1000,'k');
    set(gca,'fontweight','bold','fontsize',14);
    
    ylim([min(axis1) max(axis1)]);
    
    ylabel(feature_names{1});
    
    xlabel('Samples (thousands)');
    set(gca,'xminorgrid','on','yminorgrid','on');
    hout(1)=gca;
    plot_letter_label('a)');
    
    grid on
    poss=get(gca,'pos');
    poss(1)=0.200 ;
    set(gca,'pos',poss);
    
    %%%%%%%%%%%%
    subplot(2,2,2)
end
if scale_chc==2  %log10
    rangee=0:0.5:ceil(max(max((log10(N)))));
    if contour_chc==1
        imagesc(axis2,axis1,log10(N));
    else
        contourf(axis2,axis1,log10(N),rangee);
    end
    
else
    rangee=1:1:ceil(max(max(N)));
    if contour_chc==1
        imagesc(axis2,axis1,N);
    else
        contourf(axis2,axis1,N,rangee);
    end
    
end

hout(2)=gca;
set(gca,'fontweight','bold','fontsize',14);
set(gca,'xminorgrid','on','yminorgrid','on');

hh=colorbar('east');
%set(hh,'fontweight','bold','fontsize',14,'DecorationColor','w','DecorationColor','k');
set(hh,'fontweight','bold','fontsize',14);

poss=get(hh,'pos');
poss(1)=.94;
set(hh,'pos',poss);
cmap=interp1(linspace(0,1,(size(colormap,1))),colormap,rangee(2:end)/max(rangee));
colormap(cmap);
caxis([min(rangee) max(rangee)]);
if ~no_margin
    hhhh=plot_letter_label('b)');set(hhhh,'color','w');
end
axis('xy');
%xlabel(label2);
%ylabel(label1);
set(gca,'xminorgrid','on','yminorgrid','on');
title(sprintf('N: %i',size(X,2)));

if ~no_margin
    %%%%%%%%%%
    subplot(2,2,4)
    
    [N2,X2]=hist(X(2,:),axis2);
    
    bar(X2,N2/1000,'k');
    set(gca,'fontweight','bold','fontsize',14);
    
    hout(3)=gca;
    xlabel(feature_names{2});
    ylabel('Samples (thousands)');
    grid on
    
    xlim([min(axis2) max(axis2)]);
    plot_letter_label('c)');
    
    set(gca,'xminorgrid','on','yminorgrid','on');
    
    poss=get(gca,'pos');
    poss(2)=.2;
    set(gca,'pos',poss);
end

orient landscape

printname=sprintf('%s_vs_%s',feature_names{1},feature_names{2});


    function hhhh=plot_letter_label(charr)
        hhhh=text(0.8,0.9,charr,'units','norm','fontweight','bold','fontsize',14,'color','k');
    end
end