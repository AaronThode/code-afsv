function [N2D,printname,Ibin,hout,hprint]=hist2D(X,axis1,axis2,feature_names,contour_chc,scale_chc,trim_chc,no_margin)
%function [N,printname,Ibin]=hist2D(X,axis1,axis2,feature_names,contour_chc,scale_chc,trim_chc)
%  Plot 2D histograms of several variables...
%
%  Input:
%%% X: [Ndim Nsamples] data sample with Ndim variables (Ndim>1);
%%%  axis1, axis2: desired histogram bins.  axis2 can be a cell matrix for
%%%             Ndim> 2
%%%  feature_names: cell array of names of features for plotting
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

Nrows=size(X,1)-1;
if Nrows<1
    error('X needs to have more than one row');
end
Ibin=zeros(size(X));

if Nrows>1 & ~iscell(axis2)
    error('axis1 needs to be a cell matrix if size(X,1)>2');
elseif ~iscell(axis2)
    tmp=axis2;clear axis2
    axis2{1}=tmp;
end

if Nrows~=length(axis2)
    error('Dimensions of axis2 must match rows of X');
end

if ~exist('contour_chc')|isempty(contour_chc)
    contour_chc=1;
end
if ~exist('scale_chc')||isempty(scale_chc)
    scale_chc=2;
end

if ~exist('trim_chc')||isempty(trim_chc)
    trim_chc=1;
end

if ~exist('no_margin')||isempty(no_margin)
    no_margin=false;
end

hprint=figure;


for Irows=1:Nrows
    
    if Irows==1 && ~iscell(axis2)
        [N2D{Irows}, Ibin]=make2Dhist(axis1,axis2,X([1 Irows+1],:),trim_chc);
        
    else
        [N2D{Irows}, Ibin]=make2Dhist(axis1,axis2{Irows},X([1 Irows+1],:),trim_chc);
    end
    
    %%%%%%%%%%%%%%%%Margin along side
    if ~no_margin
        subplot(Nrows+1,2,2*Irows-1);
        [N2,X2]=hist(X(1+Irows,:),axis2{Irows});
        barh((X2),(N2)/1000,'k');
        set(gca,'fontweight','bold','fontsize',12);
        
        ylim([min(axis2{Irows}) max(axis2{Irows})]);
        ylabel(feature_names{1+Irows});
        
        xlabel('Samples (thousands)');
        set(gca,'xminorgrid','on','yminorgrid','on');
        hout(2*Irows-1)=gca;
        plot_letter_label([char(97+2*Irows-2) ')']);
        
        grid on
        poss=get(gca,'pos');
        poss(1)=0.200 ;
        set(gca,'pos',poss);
        
        %%%%%%%%%%%%
        subplot(Nrows+1,2,2*Irows);
    else
        subplot(Nrows,1,Irows)
    end
    
    
    if scale_chc==2  %log10
        rangee=0:0.5:ceil(max(max((log10(N2D{Irows})))));
        if contour_chc==1
            imagesc(axis1,axis2{Irows},log10(N2D{Irows}));
        else
            contourf(axis1,axis2{Irows},log10(N2D{Irows}),rangee);
        end
        
    else
        rangee=1:1:ceil(max(max(N2D{Irows})));
        if contour_chc==1
            imagesc(axis1,axis2{Irows},N2D{Irows});
        else
            contourf(axis1,axis2{Irows},N2D{Irows},rangee);
        end
        
    end
    colormap(jet);
    hout(2*Irows)=gca;
    set(gca,'fontweight','bold','fontsize',12);
    set(gca,'xminorgrid','on','yminorgrid','on');
    
    hh=colorbar('east');
    %set(hh,'fontweight','bold','fontsize',12,'DecorationColor','w','DecorationColor','k');
    set(hh,'fontweight','bold','fontsize',12);
    
    poss=get(hh,'pos');
    poss(1)=.94;
    set(hh,'pos',poss);
    cmap=interp1(linspace(0,1,(size(colormap,1))),colormap,rangee(2:end)/max(rangee));
    colormap(cmap);
    caxis([min(rangee) max(rangee)]);
    if Nrows>1
        %hhhh=plot_letter_label('b)');
        if ~no_margin
            hhhh=plot_letter_label([char(97+2*Irows-1) ')']);set(hhhh,'color','w');
        else
            hhhh=plot_letter_label([char(97+Irows-1) ')']);set(hhhh,'color','w');
            
        end
    end
    
    if no_margin
        ylabel(feature_names{1+Irows});
        
    end
    axis('xy');
    %xlabel(label2);
    %ylabel(label1);
    set(gca,'xminorgrid','on','yminorgrid','on');
    if Irows==1
        title(sprintf('N: %i',size(X,2)));
    end
    
    
end

if ~no_margin
    %%%%%%%%%%Reference variable
    subplot(Nrows+1,2,(2*Nrows+2))
    
    [N1,X1]=hist(X(1,:),axis1);
    bar(X1,N1/1000,'k');
    set(gca,'fontweight','bold','fontsize',12);
    
    hout(end+1)=gca;
    xlabel(feature_names{1});
    ylabel('Samples (thousands)');
    grid on
    
    xlim([min(axis1) max(axis1)]);
    %plot_letter_label('c)');
    hhhh=plot_letter_label([char(97+2*Irows) ')']);
    set(gca,'xminorgrid','on','yminorgrid','on');
    
    %poss=get(gca,'pos');
    %poss(2)=.15;
    %set(gca,'pos',poss);
end


orient landscape

printname=sprintf('%s_vs_%s',feature_names{1},feature_names{2});


    function hhhh=plot_letter_label(charr)
        hhhh=text(0.8,0.9,charr,'units','norm','fontweight','bold','fontsize',12,'color','k');
    end
end

function [N,Ibin]=make2Dhist(axis1,axis2,X,trim_chc)

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
N=N';
end %function make2Dhist