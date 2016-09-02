%%Extract a bunch of airgun samples and create a template for
%%correlation...

function airgun=extract_repetitive_template(Igun,data_all,head,param,Idebug)

airgun=[];

if isempty(Igun),
    return;
end
Icall=0;
Fs=head.Fs;

%Set up timescales
template_time_inc=0.5*3600;  %Time in hours to average airgun signals
percent_residency=0.25;
ctimes=data_all.ctime(Igun);
airgun_edge_times=min(ctimes):template_time_inc:max(ctimes);
[Nguns,airgun_bins]=histc(ctimes,airgun_edge_times);

%Set up spectrogram template
dT=param.Nfft*(1-param.ovlap)/Fs;
ImaxCall=2*round(param.energy.MaxTime/dT);
Icenter=ImaxCall/2;
Istrip=round(Fs*param.energy.bufferTime);
Tall=dT*(1:ImaxCall);

% Option 1:

for I=2:length(Nguns),

    template{I}=zeros(1+param.Nfft/2,ImaxCall);

    Iref=Igun(airgun_bins==I);
    snips_name=dir('*.snips');
    if I<length(Nguns),
        [x,nstarts_snips,npts_snips]=readEnergySnips(snips_name.name, Iref,'short','cell','keep_open');
    else
        [x,nstarts_snips]=readEnergySnips(snips_name.name, Iref,'short','cell');

    end
    
           
    for II=1:length(Iref),
        index=1:(length(x{II})-Istrip);
        load_time=data_all.ctime(Iref(II));
        temp=datevec(datenum(1970,1,1,0,0,load_time))
        fnn=dir([param.energy.file_dir '/*' num2str(temp(2)) num2str(temp(3))  '*sio']);
        fn=[param.energy.file_dir '/' fnn.name];
        [y,t,head]=readsiof(fn,load_time-1,3);
        [B0,F0,T0]=specgram(y-32768,param.Nfft,head.Fs,[],round(param.ovlap*param.Nfft));
         B0=20*log10(abs(B0+eps));



        [B,F,T]=specgram(x{II}(index),param.Nfft,head.Fs,[],round(param.ovlap*param.Nfft));
        Nt=length(T);
        B=20*log10(abs(B+eps));


        Imedian_sample=find(T<param.eq_time);
        if length(Imedian_sample)>1,
            Bmean=median(B(:,Imedian_sample).').';
        else
            Bmean=B(:,Imedian_sample);
        end
        Beq=B-Bmean*ones(1,size(B,2));
        Beq(:,Imedian_sample)=0;
        Iminn=find(Beq<0);
        Beq(Iminn)=0;
        [SEL,Imax]=max(sum(Beq));

        index2=(-floor(Nt/2):floor(Nt/2))+Icenter-(Imax-floor(Nt/2)-1);
        template{I}(:,index2)=template{I}(:,index2)+Beq;

        if Idebug>0,
            
            figure(9);
            subplot(2,1,1);
            imagesc(T,F,B);caxis([20 110]);colorbar;axis('xy');
            title(sprintf('%i of %i signals incorporated',II,length(Iref)));
            subplot(2,1,2);
            imagesc(T0,F0,B0);caxis([20 110]);colorbar;axis('xy');
            title(sprintf('Loaded directly from sio file'));
            
            figure(10);
            subplot(2,1,1);
            plot(x{II}(index));%colorbar;axis('xy');
            title(sprintf('%i of %i signals incorporated',II,length(Iref)));
            subplot(2,1,2);
            plot(t,y-32768);%colorbar;axis('xy');
            title(sprintf('Loaded directly from sio file')); 

            figure(8);
            subplot(2,1,1);
            imagesc(T,F,B);colorbar;axis('xy');
            title(sprintf('%i of %i signals incorporated',II,length(Iref)));
           
            subplot(2,1,2);
            imagesc(Tall,F,template{I});colorbar;axis('xy');
            title(sprintf('%i of %i signals incorporated',II,length(Iref)));
            pause;
        end

    end %II, end of averaging period
    template{I}=sign(template{I}/length(Iref)-percent_residency);
    if Idebug>1,
        figure(8);
        subplot(2,1,1);
        imagesc(T,F,template{I});colorbar;axis('xy');
        pause;
    end
end  %I -- calls

airgun.template=template;
airgun.edge_times=airgun_edge_times;
