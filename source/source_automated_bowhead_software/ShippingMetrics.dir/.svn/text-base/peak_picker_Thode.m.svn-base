%function peaks=peak_picker_Thode(PdB,F,df_search,frange,thresholddB,Idebug)
% Input:
% PdB: a power spectral density, in terms of dB units.  Can be a row or column vector
% F:   vector of frequencies that accompanies
% df_search:  Half-bandwidth to permit for a peak search, Hz.
%       If a vector, element I of df_search is returns in peaks{I}
% frange:  [fmin fmax] frequency range in Hz to search for peaks.
% thresholddB:  SNR (dB) of peak compared to adjacent bins
% Idebug:  Optional structure for plotting.  If exists, has following fields:
%   .fig:  Figure number.  The vector should be the same length as df_search input vector.
%   .title:  Title to paste on top
%
% Output:
% peaks: a cell structure with two fields... isi, and adp, for "ISI_tone" and "adaptive"
%       peaks{I} corresponds to element I of df_search...  The adpative does not make assumptions about signal
%       bandwidth.
%       .isi:
%           .F: frequencies at which a peak occurs
%           .P: linear power spectral density at each peak in F
%           .PdB: dB power spectral density at each peak in F
%           .PdBall: scalar with sum of linear PSD across all detected peaks (total industrial power).
%
%       .adp:
%            .F= Frequencies at which a peak occurs
%            .PdBint= integrated power across each peak bandwidth (dB)
%            pdBinttotal=total industrial power at a given time
%            .PdB= peak PSD for each peak in F
%            .bandwidth  threshold bandwidth of each peak in Hz
%


function peaks=peak_picker_Thode(PdB,F,df_search,frange,thresholddB,Idebug)

%%For debugging as a script...
%clear
%load temp
%Idebug.fig=1;
%df_search=10;
PdB=PdB(:);

P=10.^(PdB/10);  %linear units
threshold=10^(thresholddB/10);  %linear units

f_min	=	frange(1);
f_max	=	frange(2);
[~,if_min]	=	min(abs(F - f_min));
[~,if_max]	=	min(abs(f_max - F));
if_min	=	if_min(1);
if_max	=	if_max(1);
dF=F(2)-F(1);
Isearch=ceil(df_search/dF);

dP=diff(PdB(if_min:if_max));
test1=sign(dP(1:(end-1)));
test2=dP(2:end).*(dP(1:(end-1)));%The sign restricts to maxima
Ipeak=if_min+find(test1>0&test2<0);  %Note no restrictions on level


%peaks.Fmax=F(Ipeak);
%peaks.PdBmax=PdB(Ipeak);


for Ipar=1:length(Isearch)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Simple ISI_tone...%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    peaks{Ipar}.isi.F=[];
    peaks{Ipar}.isi.P=[];
    for I=1:length(Ipeak)
        index=Ipeak(I)+([-Isearch(Ipar):-1 1:Isearch(Ipar)]);
        if (P(Ipeak(I))./median(P(index))>=threshold)
            peaks{Ipar}.isi.F=[peaks{Ipar}.isi.F F(Ipeak(I))];
            peaks{Ipar}.isi.P=[peaks{Ipar}.isi.P P(Ipeak(I))];
        end
        
    end
    peaks{Ipar}.isi.PdBall=10*log10(sum(peaks{Ipar}.isi.P));
    peaks{Ipar}.isi.PdB=10*log10(peaks{Ipar}.isi.P);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%My adaptive search %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % here df_is a maximum peak width permitted
    peaks{Ipar}.adp.F=[];
    peaks{Ipar}.adp.PdBint=[];
    peaks{Ipar}.adp.PdB=[];
    peaks{Ipar}.adp.bandwidth=[];
    for I=1:length(Ipeak)
        index=Ipeak(I)+( 1:Isearch(Ipar));
        
        %%Search for upper threshold
        Iu=Ipeak(I)+find((PdB(Ipeak(I))-thresholddB)>PdB(index), 1 );
        
        index=Ipeak(I)+(-Isearch(Ipar):-1 );
        Il=Ipeak(I)-find((PdB(Ipeak(I))-thresholddB)>PdB(index), 1, 'last' );
        
        if isempty(Il)||isempty(Iu)
            continue
        end
        
        peaks{Ipar}.adp.F=[peaks{Ipar}.adp.F F(Ipeak(I))];
        peaks{Ipar}.adp.bandwidth=[peaks{Ipar}.adp.bandwidth F(Iu)-F(Il)];
        peaks{Ipar}.adp.PdBint=[peaks{Ipar}.adp.PdBint dF*sum(P(Il:Iu))];
        peaks{Ipar}.adp.PdB=[peaks{Ipar}.adp.PdB PdB(Ipeak(I))];
    end
    
    peaks{Ipar}.adp.pdBinttotal=10*log10(sum(peaks{Ipar}.adp.PdBint));
    peaks{Ipar}.adp.PdBint=10*log10(peaks{Ipar}.adp.PdBint);
    
    if exist('Idebug','var')&~isempty(Idebug)
        figure(Idebug.fig(Ipar));clf
        plot(F,PdB,'k^-',F(Ipeak),PdB(Ipeak),'ro');grid on
        
        hold on
        plot(peaks{Ipar}.isi.F,10*log10(peaks{Ipar}.isi.P),'gs');  %green square is ISI_tone estimate...
        
        
        for I=1:length(peaks{Ipar}.adp.PdB)
            Fplot=[ peaks{Ipar}.adp.F(I)+peaks{Ipar}.adp.bandwidth(I)*[-0.5 0.5]];
            Pplot=peaks{Ipar}.adp.PdB(I)*[1 1]+5;
            plot(peaks{Ipar}.adp.F,5+peaks{Ipar}.adp.PdB,'b^',Fplot,Pplot,'b-');
            
        end
        
        legend('original points','local peaks','isi_tone','adaptive')
        title(Idebug.title);
    end
end %Ipar

%end