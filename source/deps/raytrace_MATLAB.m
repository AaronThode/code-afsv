%function [x,dt,dz]=raytrace_MATLAB(ang_meas,dt_meas,ssp_chc,scenario_chc)
%  If called with no arguments, returns allowable scenarios

function [x,dt,dz]=raytrace_MATLAB(ang_meas,dt_meas,scenario_chc,ssp_chc)

x=[];
dt=[];
dz=[];

prompt = {'Water depth (m)','Source depth (m)','Profile choice (munk,linear,MURI_Feb2014)', ...
    'path length (km)','ray angles (degrees)', ...
    'Relative arrival times (msec; positive means arrives before first path)','Tilt (deg)' , ...
    'Uncertainty in timing (msec)','Uncertainty in depth (m)'};
dlg_title = 'Ray tracing parameters';
num_lines = 1;

scenarios={'MURI_Feb2014' };
%scenarios={'MURI_Feb2014', 'MURI_Feb2014_20140218T134800', 'MURI_Feb2014_20140218T134000','MURI_Feb2014_20140218T134200',  ...
    %'MURI_Feb2014_20140218T134800','MURI_Feb2014_20140218T140039', 'MURI_Feb2014_Humpback' };

if nargin==0
    x=scenarios;
    return
end

switch scenario_chc
    case 'MURI_Feb2014'
        
        def = {'4000','330','MURI_Feb2014','200',mat2str(ang_meas),mat2str(dt_meas),'-2:2','2','50'};
   
    otherwise
        def = {'4000','330','MURI_Feb2014','200',mat2str(ang_meas),mat2str(dt_meas),'-1:1','2','50'};
        %errdlg(sprintf('Scenario %s not in database!',scenario_chc));
        %return
end


answer = inputdlg(prompt,dlg_title,num_lines,def);

if isempty(answer)
    return
end
D=eval(answer{1}); %bottom depth in m
z0=eval(answer{2});
choice=answer{3};
switch choice,
    case 'munk',
        c=getcmunk(z0);
        funname='munkcalc';
        chctit='munk';
    case 'linear',
        c=getclin(z0);
        funname='linearcalc';
        chctit='linear';
    case 'MURI_Feb2014',
        c=getcgeneral_crude(z0);
        funname='generalcalc';
        chctit='MURI profile';
        
end
disp(chctit);
send=1000*eval(answer{4});
theta0=eval(answer{5});
dt_meas=eval(answer{6});
tilt=eval(answer{7});
sig_dt=eval(answer{8});  %Uncertainty in measurement, msec
sig_z=eval(answer{9});  %Uncertainty in depth

best_err=zeros(1,length(tilt));
for Itilt=1:length(tilt)
    fprintf('Tilt: %6.2f deg\n',tilt(Itilt));
    clear x r z
    theta=theta0+tilt(Itilt);
    
    Nrays=length(theta);
    theta=pi*theta/180;
    
    %Begin ray tracing
    Nsamps=0;
    for I=1:Nrays
        %Set up initial conditions
        r=0;
        z=z0;
        rho=cos(theta(I))/c;
        eta=sin(theta(I))/c;
        tau=0;
        
        x0=[r; z; rho; eta; tau];
        fprintf('Starting ray %i, launch angle %6.2f deg\n',I,theta(I)*180/pi);
        [s,x{I}]=ode45ray(funname,0,send,x0, 1e-10,0,D);
        Nsamps=max([Nsamps length(s)]);
    end
    
    clear r z
    
    dt=zeros(Nrays,Nsamps);
    dz=dt;
    cchc='krbgckrbgc';
    for I=1:Nrays
        N=size(x{I},1);
        x{I}=[x{I}; ones(Nsamps-N,1)*x{I}(end,:)];
        
        
        
        %interpolate range vs path 1
        t1=x{1}(:,5);
        r1=x{1}(:,1);
        z1=x{1}(:,2);
        
        t{I}=x{I}(:,5);
        r{I}=x{I}(:,1);
        z{I}=x{I}(:,2);
        
        Igood=find(diff(r{I})>0);
        r{I}=r{I}(Igood);t{I}=t{I}(Igood);z{I}=z{I}(Igood);
        t_intp{I}=interp1(r{I},t{I},r1);
        z_intp{I}=interp1(r{I},z{I},r1);
        
        
        %tref1=interp1(r0,t0,x{I}(:,1));
        dz(I,:)=z1'-z_intp{I}';
        dt(I,:)=t1'-t_intp{I}';  %Positive means that a given path arrives before
        %   the first path.
        
        %%%Plot data
        figure(Itilt);
        subplot(2,1,1);
        plot(x{I}(:,1)/1000,x{I}(:,2),cchc(I),'linewidth',3);
        set(gca,'fontweight','bold','fontsize',14);
        axis('ij');hold on
        xlabel('range(km)');ylabel('depth(m)');grid on
        ylim([0 D]);set(gca,'ytick',0:250:D);
        
        subplot(2,1,2);
        
        plot(r1/1000,dt(I,:)*1000,cchc(I),'linewidth',3);hold on
        set(gca,'fontweight','bold','fontsize',14);
        
        %Plot measured time
        plot([min(r1) max(r1)]/1000,dt_meas(I)*[1 1],['--' cchc(I)],'linewidth',2);
        
        xlabel('range (km)');ylabel('Time relative to first ray (msec)');
        title('Positive value means path arrives before path 1');
        
        grid on
        
        ylimm=[-150 150];
        ylim(ylimm);set(gca,'ytick',ylimm(1):10:ylimm(2));
        
        %     subplot(3,1,2);
        %     %Check raw results
        %     t0=x{1}(:,1)/1500;
        %     plot(x{I}(:,1)/1000,x{I}(:,5)-t0,cchc(I),'linewidth',3);hold on
        %     set(gca,'fontweight','bold','fontsize',14);
        %     xlabel('range (km)');ylabel('Flight time-ray1 (sec)');
        %     grid on
        
    end
    
    subplot(2,1,1);
    title(sprintf('%i rays in a %s sound speed profile, Tilt %6.2f degrees',Nrays,chctit,tilt(Itilt)));
    %title([int2str(Nrays) ' rays in a ' chctit ' sound speed profile.'])
    orient landscape
    %ylim([0 1000]);
    
    fname=sprintf('RayTimetrace_%s_tilt%3.2f.jpg',scenario_chc,tilt(Itilt));
    print(gcf,'-djpeg',fname);
    
    %save(sprintf('Ray_trace_demo_tilt_%3.2f.mat',tilt),'x','dt','dz');
    [best_err(Itilt),best_range(Itilt),best_depth(Itilt)]=error_plot(x,dt,dz,dt_meas,sig_dt,sig_z,fname,Itilt+20,1000);
    orient landscape
    [~,fname,~]=fileparts(fname);
    print(gcf,'-djpeg',sprintf('CrudeError_%s.jpg',fname));
    saveas(gcf,sprintf('CrudeError_%s.fig',fname),'fig')
end %Itilt

%%Plot error results
figure
subplot(3,1,1);
plot(tilt,best_err,'ko','linewidth',3);
set(gca,'fontweight','bold','fontsize',16);grid on;
ylabel('mean(err/uncertainty)')

subplot(3,1,2);
plot(tilt,best_range/1000,'ko','linewidth',3);
set(gca,'fontweight','bold','fontsize',16);grid on;
ylabel('Range (km)')

subplot(3,1,3);
plot(tilt,best_depth,'ko','linewidth',3);
set(gca,'fontweight','bold','fontsize',16);grid on;

xlabel('Tilt (deg)');ylabel('Depth (m)')

orient landscape
print(gcf,'-djpeg',sprintf('ErrVsTilt_RayTimeTrace_%s.jpg',scenario_chc))
end

function [best_err,rbest,zbest]=error_plot(x,dt,dz,dt_meas,sig_dt,sig_z,fname,Ifig,D)

if size(dt_meas,2)>1
    dt_meas=dt_meas';
end

Nrays=size(dt,1);
Npts=size(dt,2);
err_dt=NaN*ones(1,Npts);
err_z=err_dt;
for I=1:Npts
    err_dt(I)=sum((dt_meas-dt(:,I)*1000).^2)./sig_dt.^2;
    err_z(I)=sum((dz(:,I)).^2)./sig_z.^2;
end

err_dt=sqrt(err_dt);
err_z=sqrt(err_z);
err_all=err_dt+err_z;

[best_err,Ir]=min(err_dt);
fprintf('modeled time delays for timing error: %s\n',mat2str(dt(:,Ir)*1000));

[best_err,Ir]=min(err_all);
fprintf('modeled time delays for total error: %s\n',mat2str(dt(:,Ir)*1000));


cchc='krbgckrbgc';
figure(Ifig);
subplot(3,1,1);
for I=1:Nrays
    plot(x{I}(:,1)/1000,x{I}(:,2),cchc(I),'linewidth',3);axis('ij');hold on
   
end
set(gca,'fontweight','bold','fontsize',18);
 xlabel('range(km)');ylabel('depth(m)');grid on;ylim([0 D])
    set(gca,'ytick',0:250:D);
title(fname,'interp','none');

subplot(3,1,2);
semilogy(x{1}(:,1)/1000,(err_dt),x{1}(:,1)/1000,(err_z),'linewidth',3);
set(gca,'fontweight','bold','fontsize',18);

legend('Timing error','depth error');ylim([0.1 100]);
xlabel('range(km)');ylabel('rms err ');grid on

subplot(3,1,3);
semilogy(x{1}(:,1)/1000,(err_all),'linewidth',3);grid on;ylim([1 100]);
set(gca,'fontweight','bold','fontsize',18);
xlabel('range(km)');ylabel('total err ');grid on

rbest=x{1}(Ir,1);
zbest=x{1}(Ir,2);

fprintf('Error: %6.2g, Best range: %6.2f km, Best depth: %6.2f m\n**\n',best_err,rbest/1000,zbest);


end
