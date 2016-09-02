%function [match_index,dt,best_mismatch]=match_stations(station_anchor,station_other,params,debug_params,fname_anchor,fname_other);
%
% Input:
%   params.names={'ctime','fmin','Centroid','BoundingBox'};
%   params.weight=[0 1 0 0];
%   params.tol=[0 25 0 0 ];
%   params.time_tol=distance/1500;  %Time separation  in seconds between stations, used
%       to set a time delay cutoff.
%   params.Fs: sampling frequency in Hz.
%   station_anchor..  Reference anchor station...
%   station_other,    Anchor station to search for matches...
%   station(Istation).SNR=[1 x Ncall] vector of PSD-based SNR
%   station(Istation).ctime=best_ctimes;
%   debug_params:
%       J_anchor_start:  indicies of call in anchor station to be checked
%       cross_channel:  if 1, some debug printing, if 2, more debug
%               printing..
%   fname_anchor,fname_other:  complete filenames of raw acoustic data
%   associated with anchor and other station.
%
% Output:
%    match_index: [2 by Nmatch] matrix of linked stations, row 1 are
%           indicies of anchor station
%    Igood_other: [1 x Ncall] vector of indicies of station_other that
%           match station_anchor indicies in Igood_anchor.  If no match value of -1.
%    Igood_anchor:  indicies of original station_anchor object that survive
%       the match_station.  Note that 'survive' means that a lower SNR match
%       has been found on another station, OR no match has been found at all.
%    dt: vertical vector of arrival_time(anchor)-arrival_time(other).
%    best_mismatch:  feature error of match at Icall...


function [match_index,dt,dt_xcorr,best_mismatch]=match_stations(station_anchor,station_other,params,debug_params,fname_anchor,fname_other)

dt=[];
dNmax=0;
%anchor_feature_str='SEL';  %Feature used to establish 'anchor' station: SEL, SNR, and TotalDuration are possible choices

%%Revise feature names if matched filter option is desired...

names{1}='xcorr';
params.names=[];
params.names{1}='xcorr';
params.tol(1)=params.xcorr.tol;
params.weight(1)=1;
compare_function=str2func(params.compare_function_str);

%I_anchor_start=[37    51    55    66    93   136 137 156 239];
I_anchor_start=1;
if exist('debug_params')&&debug_params.cross_channel>1,
    % I_anchor_start=min(debug_params.I_anchor_start);
    Ncalls=length(debug_params.I_anchor_start);
    %Igood_anchor=J_anchor_start:Ncalls;
else
    Ncalls=length(station_anchor.ctime_min);
    Ncalls=200;
end

%%%%Initialize...
match_index=zeros(2,Ncalls);
match_index(1,:)=1:Ncalls;
dt=NaN*match_index(1,:);
dt_xcorr=dt;
best_mismatch=dt;
Icount=0;

scratchpad=spalloc(Ncalls+5,5+length(station_other.ctime_min),5000);
%%March through each call on the anchor station
%for Icall=I_anchor_start:Ncalls
for Icall=I_anchor_start:Ncalls
    %Icall=match_index(I,1);
    if rem(Icall, 200)==0, disp(sprintf('%6.2f percent complete...',round(100*Icall/Ncalls)));end
    best_mismatch(Icall)=Inf;  %lowest mismatch for each "other" shape number

    anchor_time=station_anchor.ctime_min(Icall);

    %%%
    %Are detections within allowable time delay?
    %%%%
    dt_test=(anchor_time-station_other.ctime_min);
    %Idt_test=find(abs(dt_test)<=params.time_tol);
    %Icandidate=Irange(Idt_test);
    Icandidate=find(abs(dt_test)<=params.time_tol);
    if isempty(Icandidate)
        %Irange=Irange+1;
        if debug_params.cross_channel>0,
            disp(sprintf('No signals pass %10.6f sec time tolerance',params.time_tol));
        end
        continue;
    end

    %Compare feature list
    mismatch_err=zeros(1,length(Icandidate));
    if debug_params.cross_channel==2,
        disp('Cycling through combinations');
        bufferTime=0;
        debug_spectrogram_display(station_anchor,station_other,Icall,Icandidate,fname_anchor,fname_other,bufferTime,debug_params);
    end


    %%Select best-fit other match...
    tol_test=Inf*ones(1,length(Icandidate));
    for L=1:length(Icandidate)
        if (scratchpad(Icall,Icandidate(L))==0)
            [tol_test(L),dN(L)]=image_match(compare_function,station_anchor,station_other,Icall,Icandidate(L),debug_params.cross_channel);
            scratchpad(Icall,Icandidate(L))=tol_test(L);
        else
            tol_test(L)=scratchpad(Icall,Icandidate(L));
            dN(L)=0;
        end
    end
    Igood=find(tol_test<=params.tol);
    Icandidate=Icandidate(Igood);
    tol_test=tol_test(Igood);
    %Idt_test=Idt_test(Igood);
    dN=dN(Igood);

    [tol_test,Isort]=sort(tol_test,'ascend');  %Lower score is best and comes first
    Icandidate=Icandidate(Isort);
    %Idt_test=Idt_test(Isort);
    dN=dN(Isort);

    if debug_params.cross_channel==2,
        print_debug_summary;
    end

    %%If no matches, keep this as an anchor--other stations may have it.
    if isempty(Icandidate)
        if debug_params.cross_channel==2,
            disp(sprintf('No matches pass tolerance criteria.\n\n\n\'));
        end

        continue
    end

    %Reciprical test

    if scratchpad(Icall)==Icandidate(1)  %If a previous check yielded this index..
        Ibest=1;
    else
        Ibest=[];
        minval_min=Inf;
        % tol_test2=Inf*ones(1,length(Icandidate));
        for K=1:length(Icandidate)
            dt_test2=(station_other.ctime_min(Icandidate(K))-station_anchor.ctime_min);
            Icandidate2=find(abs(dt_test2)<=params.time_tol);
            %Icandidate2=Icandidate2(Icandidate2>=Icall);  %Need to permit
            %       a superior match to endure...
            tol_test2=Inf*ones(1,length(Icandidate2));
            for L=1:length(Icandidate2)
                if Icandidate2(L)==Icall
                    tol_test2(L)=tol_test(K);
                elseif (scratchpad(Icandidate2(L),Icandidate(K))==0)
                    [tol_test2(L)]=image_match(compare_function,station_other,station_anchor,Icandidate(K),Icandidate2(L),debug_params.cross_channel);
                    scratchpad(Icandidate2(L),Icandidate(K))=tol_test2(L);
                else
                    tol_test2(L)=scratchpad(Icandidate2(L),Icandidate(K));
                end

            end
            %[tol_test2]=image_match(compare_function,station_other,station_anchor,Icandidate(K),Icandidate2,debug_params.cross_channel);
            [minval,Ival]=min(tol_test2);
            if minval<=min([params.tol minval_min])&Icandidate2(Ival)==Icall
                Ibest=K;
                minval_min=minval;

            end


        end
    end

    if isempty(Ibest)
        %disp('No countermatch found');
        continue
    else
        %disp('Match found');

    end

    %[low_mismatch,Ibest]=min(mismatch_err);


    %%If a match, get 'anchor eature of matching station
    %anchor_feature_other=station_other.(anchor_feature_str)(Icandidate(Ibest));
    ctime_match=station_other.ctime_min(Icandidate(Ibest));


    %%Assign the best_match information...
    match_index(2,Icall)=Icandidate(Ibest);
    best_mismatch(Icall)=tol_test(Ibest);
    dt_xcorr(Icall)=dt_test((Ibest))+params.dT*dN(Ibest);
    dt(Icall)=dt_test((Ibest));
    dNmax=max([ abs(dN(Ibest)*params.dT) params.dT*dNmax]);
    %dt(Icall)=dt_test(Icandidate(Ibest));

    %hold the horses if we have a repeat...
    testt=find(Icandidate(Ibest)==match_index(2,:));
    if length(testt)>1
        disp(sprintf('other station index %i has been selected by Icalls %s with best mismatches %s', ...
            Icandidate(Ibest),mat2str(match_index(1,testt)),mat2str(best_mismatch(1,testt))));

    end

    if debug_params.cross_channel>0,
        print_debug_match_summary;
    else
        Icount=Icount+1;
        % disp(sprintf('%ith match made on call %i',Icount,Icall));
    end


end %Icall


disp(sprintf('Maximum dt from spectrogram adjustment is %10.8f',dNmax));
%%Clean up results

Ipass=find(match_index(2,:)>0);
match_index=match_index(:,Ipass);
dt=dt(Ipass);
dt_xcorr=dt_xcorr(Ipass);
best_mismatch=best_mismatch(Ipass);

%%%%%%%%%Inner functions%%%%%%%%%%%%%%%%%%%%

    function print_debug_match_summary
        disp(sprintf('New best match for Call %i found on other station at index %i: ',Icall,Igood_other(Icall)));
        disp(sprintf('Tolerance criteria: %10.6f\n',params.tol'));

        %         for If=1:length(params.names)
        %             switch params.compare_function_str
        %                 case 'image_match'
        %                     tol_test=image_match(station_anchor,station_other,Icall,Igood_other(Icall),0);
        %                     Igood=find(tol_test<=params.tol(If));
        %                     disp(sprintf(' Feature %10s:\t weight %6.2f: anchor value: %6.2f, best match value: %6.2f err: %6.4f', ...
        %                         names{If},params.weight(If),anchor_feature,test_feature,tol_test));
        %                 case 'matched_filter'  %we are to cross-correlate image.
        %                     [tol_test,test_feature,anchor_feature]=matched_filter(params.Fs,station_anchor,station_other,Icall,Igood_other(Icall),fname_anchor,fname_other,bufferTime,debug_params);
        %                     Igood=find(tol_test<=params.tol(If));
        %                     disp(sprintf(' Feature %10s:\t weight %6.2f: anchor value: %6.2f, best match value: %6.2f err: %6.4f', ...
        %                         names{If},params.weight(If),anchor_feature,test_feature,tol_test));
        %                 case 'spectrogram_correlation'  %we are to cross-correlate image.
        %                     anchor_feature=double(logical(full(station_anchor.feature.Image{Icall})));
        %                     test_feature=double(logical(full(station_other.feature.Image{Igood_other(Icall)})));
        %                     tol_test=image_correlate(anchor_feature,test_feature,debug_params.cross_channel);
        %                     disp(sprintf(' Feature %10s:\t weight %6.2f:  err: %6.4f', ...
        %                         params.names{If},params.weight(If),tol_test));
        %                 otherwise
        %                     anchor_feature=station_anchor.feature.(names{If})(Icall);
        %                     test_feature=station_other.feature.(names{If})(Igood_other(Icall));
        %                     tol_test=compare_function(test_feature,anchor_feature);
        %                     %tol_test=abs(test_feature-anchor_feature)./anchor_feature;
        %                     disp(sprintf(' Feature %10s:\t weight %6.2f: anchor value: %6.2f, best match value: %6.2f err: %6.4f', ...
        %                         params.names{If},params.weight(If),anchor_feature,test_feature,tol_test));
        %
        %             end
        %
        %         end
        disp(sprintf('dt: %10.6f',dt(Icall)));
        disp(sprintf(' Local time: %s',ctime2str(ctime_match)));
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

        %locations_all{indicies}=locations;

        disp(' ');
        %pause

    end

%%%%%%%print_debug_summary
    function print_debug_summary
        disp(' ');
        disp('%%%%%%%%%%%match_stations.m%%%%%%%%%%%%%%%%%%');
        disp(sprintf('Matching call %i:',Icall));
        disp(sprintf('Feature %10s:\t Weight: %6.2f, Tolerance criteria: %10.6f',params.names{If},params.weight(If),params.tol(If)));

        disp(sprintf('station index: %i anchor values: %6.2f feature: %6.2f err: %10.4f dt: %10.6f SNR: %10.6f\n', ...
            [Icandidate; anchor_feature*ones(size(Icandidate)); test_feature*ones(size(Icandidate)); tol_test; dt_test(Icandidate); station_other.SNR(Icandidate);]));

        if ~isempty(Igood),
            possible_ctime=station_other.ctime_min(Icandidate(Igood));
            for KK=1:length(Igood),
                disp(sprintf('Candidate(s) %i , other station index %i at %s pass this round.',Igood(KK)',Icandidate(Igood(KK))', ...
                    datestr(datenum(1970,1,1,0,0,possible_ctime(KK))',31)));
            end
        end
        disp(' ');

        pause

    end
end




%%%%%%%%image_match%%%%%%%%%%%%
function [tol_test,dN_cand]=image_match(similarity_function,station_anchor,station_other,Icall,Icandidate,Iplot)
%funcstr=func2str(similarity_function);
tol_test=[];
if isempty(Icall),
    return;
end

if isempty(Icandidate)
    return;
end

anchor_image=full(station_anchor.Image{Icall});

Nsegs=station_anchor.Nsegments(Icall);

tol_test=Inf*ones(Nsegs,length(Icandidate));
dN_cand=tol_test;

%%Images may be segmented...
for I=1:Nsegs,
    anchor_subimage=(ismember(anchor_image,I));
    for J=1:length(Icandidate),
        other_image=full(station_other.Image{Icandidate(J)});
        Nsegs_other=station_other.Nsegments(Icandidate(J));
        dist=zeros(1,Nsegs_other);
        dist_min=dist;
        dN=dist;
        for K=1:Nsegs_other
            other_subimage=full(ismember(other_image,K));
            [dist(K),dN(K),dist_min(K),max_ovlap]=correlate_image(anchor_subimage, other_subimage,Iplot);

            %tol_test(1,J)=min([tol_test(1,J) dist(K)]);

            if 1==1;
                figure(1)
                subplot(2,2,1)
                imagesc(anchor_image);
                title(int2str(Icall));
                subplot(2,2,2)
                imagesc(other_image);
                title(int2str(Icandidate(1)));
                subplot(2,2,3)
                imagesc(anchor_subimage);
                subplot(2,2,4);
                imagesc(other_subimage);
                title(sprintf('log corr: %8.4f, standard corr: %8.4f max_ovlap %i',dist(K),dist_min(K),round(max_ovlap)));
                pause

            end
        end
        [tol_test(I,J),Imin]=min(dist);
        dN_cand(I,J)=dN(Imin);

    end


end

if Nsegs>1
    disp('WARNING! NOT DEBUGGED!');
    keyboard;
    [tol_test,Imin]=min(tol_test);
    dN_cand=dN_cand(Imin,:);

end

end

%%%%%%%%%debug_spectrogram_display.m%%%%%%%%%
function debug_spectrogram_display(station_anchor,station_other,Icall,Icandidate,fname_anchor,fname_other,bufferTime,debug_params)

if isempty(Icall),
    return;
end

if isempty(Icandidate)
    return;
end
figure(1);%set(gcf,'pos',[7   228   590   832]);
[y,start_ctime]=extract_signal_from_station(station_anchor,Icall,fname_anchor,'debug');
debug_params.morph.equalization(:,1)=station_anchor.param.equalization_freq;
debug_params.morph.equalization(:,2)=station_anchor.equalization(:,Icall);
figure(1);
[SS,FF,TT,PP]=spectrogram(y,hanning(debug_params.Nfft),round(debug_params.ovlap*debug_params.Nfft),debug_params.Nfft,debug_params.Fs,'yaxis');
subplot(2,3,1)
imagesc(TT,FF,10*log10(PP));axis('xy');
xlimm=xlim;
xlimm(1)=debug_params.morph.eq_time;
xlim(xlimm);
caxis([40 120]);
set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');

subplot(2,3,2);
[stats,final_image,Beq,raw_image]=extract_image_features(y,start_ctime,debug_params,2);
if ~isempty(stats)
    FF=stats(1).dF*(0:(size(final_image,1)-1));
    TT=stats(1).dT*(0:(size(final_image,2)-1));
    plot_morph_image(TT,FF,Icall,[],final_image);
else
    disp('Funny, this time around couldn''t get a morph result');
    keyboard;
end
title(sprintf('%s start: %s %ith Shape at %6.2f sec',fname_anchor,ctime2str(start_ctime),Icall,station_anchor.ctime_min(Icall)-start_ctime-debug_params.morph.eq_time));

subplot(2,3,3)
imagesc(TT,FF,station_anchor.Image{Icall});
title('stored image');

for I=1:length(Icandidate)
    [y,start_ctime]=extract_signal_from_station(station_other,Icandidate(I),fname_other,'debug');
    [SS,FF,TT,PP]=spectrogram(y,hanning(debug_params.Nfft),round(debug_params.ovlap*debug_params.Nfft),debug_params.Nfft,debug_params.Fs,'yaxis');


    figure(1);
    subplot(2,3,4)
    imagesc(TT,FF,10*log10(PP));axis('xy');
    xlimm=xlim;
    xlimm(1)=debug_params.morph.eq_time;
    xlim(xlimm);
    set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
    caxis([40 120]);

    subplot(2,3,5);
    debug_params.morph.equalization(:,1)=station_other.param.equalization_freq;
    debug_params.morph.equalization(:,2)=station_other.equalization(:,Icandidate(I));

    [stats,final_image2,Beq,raw_image]=extract_image_features(y,start_ctime,debug_params,2);

    if ~isempty(stats)
        FF=stats(1).dF*(0:(size(final_image2,1)-1));
        TT=stats(1).dT*(0:(size(final_image2,2)-1));
        plot_morph_image(TT,FF,Icall,[],final_image2);
    else
        disp('Funny, this time around couldn''t get a morph result');
        keyboard;
    end
    title(sprintf('%s, Candidate %i of %i, start: %s Shape at %6.2f sec',fname_other,Icall, length(Icandidate), ...
        ctime2str(start_ctime),station_other.ctime_min(Icandidate(I))-start_ctime-debug_params.morph.eq_time));

    subplot(2,3,6)
    imagesc(TT,FF,station_other.Image{Icandidate(I)});
    pause;

    %[dist]=similarity_function(final_image, final_image2);
    %disp(sprintf('Best hausdorf distance is %6.2f',min(dist{1})));
    %keyboard;

end


end


%%%%%%plot_morph_image.m%%%%%%%%%%
%%%plot morphological image output

function plot_morph_image(TT,FF,station_index,index_offset_t,my_image)

imagesc(TT,fliplr(FF),sum(my_image,3))
axis('xy')
set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');

for Index=1:length(station_index),
    try,
        line(index_offset_t(Index)*[1 1],[0 max(FF)],'color',[1 1 1]);
        text(index_offset_t(Index),10,num2str(station_index(Index)),'color',[1 1 1]);
    catch
    end
end
end

function [HHlog,dN,HH,max_ovlap]=correlate_image(Y, X,Iplot)

HH=Inf;dN=0;HHlog=Inf;max_ovlap=[];

X=double(X); Y=double(Y);
M=max([size(X,2) size(Y,2)]);
Igood_x=find(sum(X,2)>0);
Igood_y=find(sum(Y,2)>0);
Igood=intersect(Igood_x, Igood_y);
if isempty(Igood);
    return
end


% 
% XX=[X(Igood,:) zeros(length(Igood),M-size(X,2))]';
% YY=[Y(Igood,:) zeros(length(Igood),M-size(Y,2))]';
% Xmag=sum(sum(X));
% Ymag=sum(sum(Y));
% max_pixels=min([Xmag Ymag]);
% 
% Irevx=ones(size(XX,1),1)*((-1).^(1:length(Igood)));
% %Irevy=((-1).^(1:length(Igood))')*ones(1,size(Y,2));
% XX=XX.*Irevx;YY=YY.*Irevx;
% %XX=(X(Igood,:).*Irevx)';YY=(Y(Igood,:).*Irevy)';
% [tmp1,lags]=xcorr(XX(:),YY(:),M);
% [max_ovlap,Ioffset]=max(tmp1);
% max_ovlap=max([0 max_ovlap]);
% dN=lags(Ioffset);
% 
% HH=1-(max_ovlap/max_pixels);

%%Add LOG function over binary image

XX1=Logg(X);
YY1=Logg(Y);
Xmag=(sum(sum((XX1.^2))));
Ymag=(sum(sum((YY1.^2))));
max_pixels2=min([Xmag Ymag]);

XX=[XX1(Igood,:) zeros(length(Igood),M-size(X,2))]';
YY=[YY1(Igood,:) zeros(length(Igood),M-size(Y,2))]';
 Irevx=ones(size(XX,1),1)*((-1).^(1:length(Igood)));

XX=XX.*Irevx;YY=YY.*Irevx;
%XX=(X(Igood,:).*Irevx)';YY=(Y(Igood,:).*Irevy)';
[tmp2,lags]=xcorr(XX(:),YY(:),M);
[max_ovlap,Ioffset]=max(tmp2);
max_ovlap=max([0 max_ovlap]);
dN=lags(Ioffset);


HHlog=1-(max_ovlap/max_pixels2);

if 1==1
    %
    XX=[X(Igood,:) zeros(length(Igood),M-size(X,2))]';
    YY=[Y(Igood,:) zeros(length(Igood),M-size(Y,2))]';
    Xmag=sum(sum(X));
    Ymag=sum(sum(Y));
    max_pixels=min([Xmag Ymag]);

    Irevx=ones(size(XX,1),1)*((-1).^(1:length(Igood)));
    %Irevy=((-1).^(1:length(Igood))')*ones(1,size(Y,2));
    XX=XX.*Irevx;YY=YY.*Irevx;
    %XX=(X(Igood,:).*Irevx)';YY=(Y(Igood,:).*Irevy)';
    [tmp1,lags]=xcorr(XX(:),YY(:),M);
    [max_ovlap,Ioffset]=max(tmp1);
    max_ovlap=max([0 max_ovlap]);
    dN=lags(Ioffset);

    HH=1-(max_ovlap/max_pixels);
    figure;
    subplot(2,2,1);
    imagesc(XX1);colorbar
    subplot(2,2,2);
    imagesc(YY1);colorbar
    subplot(2,1,2);
    plot(lags,1-tmp2./max_pixels2,'k',lags,1-tmp1./max_pixels,'r');ylim([0 1.2]);
    legend('log','standard')
    pause;close
end



end

function Xout=Logg(X)
Xout=X;
Xsum=sum(X);
Icol=find(Xsum>0);
for I=Icol
    Igood=find(X(:,I)>0);
    dIgood=diff(Igood);
    if isempty(Igood)
        continue
    end

    Ijump=find(dIgood>2);
    Irange=sort([Igood(1) Igood(Ijump)' Igood(Ijump+1)' Igood(end)]);
    Irange=reshape(Irange,2,length(Irange)/2);
    
    for J=1:size(Irange,2)
         Imed=floor(mean(Irange(:,J)));
        BW=floor(0.5*diff(Irange(:,J)));
       
        index=((Irange(1,J)-4*BW):(Irange(2,J)+4*BW));

        x2sig2=((index-Imed)./BW).^2;

        Xout(index,I)=(1-x2sig2).*exp(-0.5*x2sig2);
    end
end


%Xout=imfilter(X,h,'replicate');
%keyboard;
end

% function Xout=Logg(X)
% Xout=X;
% Xsum=sum(X);
% Icol=find(Xsum>0);
% for I=Icol
%     Igood=find(X(:,I)>0);
%     dIgood=diff(Igood);
%     if isempty(Igood)
%         continue
%     end
% 
%     Ijump=find(dIgood>2);
%     Irange=sort([Igood(1) Igood(Ijump)' Igood(Ijump+1)' Igood(end)]);
%     Irange=reshape(Irange,2,length(Irange)/2);
%     
%     for J=1:size(Irange,2)
%         index=(Irange(1,J):Irange(2,J));
% 
%         Imed=floor(median(index));
%         BW=floor(length(index)/2);
%         x2sig2=((index-Imed)./BW).^2;
%         Xout(index,I)=(1-x2sig2).*exp(-0.5*x2sig2);
%     end
% end
% 
% 
% %Xout=imfilter(X,h,'replicate');
% %keyboard;
% end

