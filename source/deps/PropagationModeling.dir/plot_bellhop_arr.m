function [Arr,Pos]=plot_bellhop_arr(fname)
%function [arr,pos]=plot_bellhop_arr(fname)
%  Plot arrival times, arrival angles, and amplitudes of a bellhop asc ARR
%  file.
close all;

dt=0.0002; %0.2 msec spacing;

Narrmx = 10;
if isempty(strfind(fname,'.arr'))
    fname=[fname '.arr'];
end

[ Arr, Pos ] = read_arrivals_asc( fname, Narrmx );


z_source=Pos.s.depth;
z_rcvr=Pos.r.depth;
r_rcvr=Pos.r.range;

nr_z=length(z_rcvr);
nr_r=length(r_rcvr);

pchc=[];
for I=1:20
    pchc=[pchc 'bgrcmyk'];
end
%lchc=['-------:::::::'];
max_plot=min([length(r_rcvr) 4]);
Iplot=0;
figure

for Ir=1:nr_r  %For every receiver range, sort times
    
    Iplot=Iplot+1;
    if Iplot>max_plot,figure,Iplot=1;end
    subplot(max_plot,3,1+3*(Iplot-1));
    
    %%% Construct time axis...
    time=squeeze(Arr.delay(Ir,:,:));
    time_reduced=NaN*ones(size(time));
    amp=-20*log10(abs(squeeze(Arr.A(Ir,:,:))));
    names=fieldnames(Arr);
    names=names(1:6);
    
    
    for Iz=1:nr_z
        
        Igood= ~isinf(amp(:,Iz));
        time_reduced(Igood,Iz)=(time(Igood,Iz)-min(time(Igood,Iz)));  %msec
        [time_reduced(:,Iz),Isort]=sort(time_reduced(:,Iz));
        amp(:,Iz)=amp(Isort,Iz);
        for In=1:length(names)
            Arr_sort.(names{In})(Ir,:,Iz)= Arr.(names{In})(Ir,Isort,Iz);
        end
        
        
        %plot(time,amp,['o' pchc(Ir)]);
        %line([1;1]*time,[zeros(1,length(Igood));amp],'color',pchc(Ir),'linestyle',lchc(Ir));
    end
    
    tarr=0:dt:(dt+max(max(time_reduced)));
    TA=NaN*ones(nr_z,length(tarr));
    It= 1+round(time_reduced./dt);
    for Iz=1:nr_z
        Igood= ~isinf(amp(:,Iz));
        
        TA(Iz,It(Igood,Iz))=amp(Igood,Iz);
        
    end
    
    if min(size(TA)>1)
        h = pcolor(tarr*1000,z_rcvr,TA );shading flat
    else
        plot(tarr*1000,TA);grid on
    end
    grid on
    xlimm=xlim;
    title(sprintf('r_rcvr = %6.2f m',r_rcvr(Ir)));
    xlabel('msec');ylabel('Source depth (mm)');
    colorbar('southoutside')
    
    
    subplot(max_plot,3,2+3*(Iplot-1));
    for Iz=1:nr_z
        Igood= ~isinf(amp(:,Iz));
        
        for Iray=1:sum(Igood)
            switch Arr_sort.NumTopBnc(Ir,Iray,Iz)
                case 0
                    sym='o';
                case 1
                    sym='x';
                case 2
                    sym='s';
                otherwise
                    sym='o';
            end
            plot(1000*time_reduced(Iray,Iz),Arr_sort.SrcAngle(Ir,Iray,Iz),[sym pchc(Iz)]);
           
            hold on
        end
        grid on
        
    end
    axis('ij');
    xlabel('Reduced time (msec)');
    ylabel('SrcAngle (deg)');
    %legend(num2str(r_rcvr,6));
    xlim(xlimm);
    
    subplot(max_plot,3,3+3*(Iplot-1));
    for Iz=1:nr_z
        Igood= ~isinf(amp(:,Iz));
        
        for Iray=1:sum(Igood)
            switch Arr_sort.NumTopBnc(Ir,Iray,Iz)
                case 0
                    sym='o';
                case 1
                    sym='x';
                case 2
                    sym='s';
                otherwise
                    sym='o';
            end
            plot(1000*time_reduced(Iray,Iz),Arr_sort.RcvrAngle(Ir,Iray,Iz),[sym pchc(Iz)]);
            hold on
        end
        grid on
        
    end
    axis('ij')
    xlabel('Reduced time (msec)');
    ylabel('RcvrAngle (deg)');
    %legend(num2str(r_rcvr,6));
    xlim(xlimm);
    
    
    
end %Ir

%legend(num2str(r_rcvr,6));hold on

end