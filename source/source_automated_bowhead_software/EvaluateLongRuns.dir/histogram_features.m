%function histogram_features(auto,auto_ctime, Iref,param, goodFile);

function [Iref_out,Ipass]=histogram_features(auto,auto_ctime, Iref, goodFile);
%load temp

Iref_out=[];
Ipass=[];
yes=menu('Feature stripping?','Yes','No');

Iref_out=Iref;
if yes==2
    return
end

redo=3;
while redo==3
    featurenames=fieldnames(auto.locations{1}{1}.feature);
    Nfeature0=length(featurenames);
    featurenames{(Nfeature0+1)}='Totalduration'; Nfeature=Nfeature0+1;
    featurenames{(Nfeature+1)}='Totalfmax'; Nfeature=Nfeature+1;
    featurenames{(Nfeature+1)}='Totalfmin'; Nfeature=Nfeature+1;
    featurenames{(Nfeature+1)}='Contour_Area'; Nfeature=Nfeature+1;
    featurenames{(Nfeature+1)}='SEL'; Nfeature=Nfeature+1;
    featurenames{(Nfeature+1)}='SNR'; Nfeature=Nfeature+1;
    
    Ichc(1)=menu('Select first feature:',featurenames);
    Ichc(2)=menu('Select second feature:',featurenames);
    
    prompt={'Number of histogram bins:'};
    name='Bin number';
    numlines=1;
    defaultanswer={'100'};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    Nbins=str2num(answer{1});
    
    features=zeros(2,2*length(Iref));
    Ifeature=zeros(1,2*length(Iref));
    myfeaturenames=featurenames(Ichc);
    Icol=1;
    for I=1:length(Iref)
        Iplot_ref=Iref(I);
        Istations=(auto.locations{1}{Iplot_ref}.station_indicies>0);
        
        for J=1:2
            if Ichc(J)<=Nfeature0
                tmp=[auto.locations{1}{Iplot_ref}.feature(Istations).(myfeaturenames{J})];
                
            else
                tmp=[auto.locations{1}{Iplot_ref}.(myfeaturenames{J})(Istations)]';
            end
            
            features(J,Icol:(Icol+length(tmp)-1))=tmp;
            Ifeature(Icol:(Icol+length(tmp)-1))=I;
            
        end %J
        Icol=Icol+length(tmp);
        
    end
    
    for J=1:2
        axiss{J}=linspace(min(features(J,:)),max(features(J,:)),Nbins);
    end
    
    [N,printname,Ibin]=hist2D(features,axiss{1},axiss{2},myfeaturenames);
    
    disp('Select two points to "filter" space, or just hit return to skip:');
    tmp=ginput(2);
    if isempty(tmp)
        return
    end
    
    xlimm=sort(tmp(:,1));ylimm=sort(tmp(:,2));
    
    %Pass localizations that have at least one event meeting feature criteria...
    fprintf('Boundaries x: %6.2f, Boundaries y: %6.2f\n',xlimm,ylimm);
    Ipass=find(features(1,:)>=ylimm(1)&features(1,:)<=ylimm(2)&features(2,:)>=xlimm(1)&features(2,:)<=xlimm(2));
    Ipass=unique(Ifeature(Ipass));
    
    subplot(2,2,2);
    rectangle('position',[xlimm(1) ylimm(1) diff(xlimm) diff(ylimm)]);
    
    redo=menu('Choice?','Retain inside box','Retain outside of box','Redo');
    
    if redo==1
        Iref_out=Iref(Ipass);
    elseif redo==2  %Return values outside of box...
        Ipass=setdiff(1:length(Iref),Ipass);
        Iref_out=Iref(Ipass);
    end
end
    
