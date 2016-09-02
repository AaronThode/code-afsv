function call=connect_segments(connect,detect,cstart,Idebug),
%% function call=connect_segments(connect,detect,cstart),
%% Connect segments into contours
%% Assume segment 1 occurs earlier than segment 2, etc...
%%
%% connect.seg_gap_t: parameter for maximum time jump allowed
%% connect.seg_gap_f: parameter for maximum frequency jump allowed.
%% connect.seg_total_duration.  minimum length of entire contour
%% detect.T{Iseg}, detect.F{Iseg}: vector of times (sec) of segment# Iseg,
%%      frequencies etc. No other parameters needed.
%% cstart: start time in c-time that relative times in T are referenced to.
%% Idebug:  If greater than one, display criteria for failure

%detect_out=detect;
Idebug=0;
%keyboard;
%Ntduration=floor(connect.seg_total_duration./(T(2)-T(1)));

Icall=0;

if length(detect.T)>1,  %If more than one segment
     %disp(sprintf('%i segments to test',length(detect.T)));
        
    for I=1:length(detect.T),
         
        if any(detect.T{I}<0), %if we have already joined this segment
            continue
        end

        [max_t,Imax]=max(detect.T{I});
        max_f=detect.F{I}(Imax);  %Frequency at end of segment
       
        temp.T=detect.T{I};
        temp.F=detect.F{I};
        %temp.score=detect.score(I);
        
        for J=(I+1):length(detect.T),

            [min_t,Imin]=min(detect.T{J});
            min_f=detect.F{J}(Imin);

             if Idebug==1,
                figure(5);
                plot(temp.T,temp.F,'kx',detect.T{J},detect.F{J},'rx');
                %keyboard;

            end
            if ((min_t-max_t)<connect.seg_gap_t)&(min_t-max_t)>0&(abs(min_f-max_f)<connect.seg_gap_f), %segments connected
                %&(temp.score+detect.score(J)>Ntduration)  %this test
                %removed until end...
                [temp.T,Isort]=sort([temp.T detect.T{J}]);
                temp.F=[temp.F detect.F{J}];
                temp.F=temp.F(Isort);
                %temp.score=temp.score+detect.score(J);

                %Reset maximum value of segment...
                [max_t,Imax]=max(detect.T{J});
                max_f=detect.F{J}(Imax);
                
                detect.T{J}=-100;
                detect.F{J}=-100;
               
            end

        end  %J
        
        %Evaluate final call--is it long enough?
       % if temp.score>Ntduration,
         if max(temp.T)-min(temp.T)>=connect.seg_total_duration&max(temp.T)-min(temp.T)<=connect.seg_max_duration,
            Icall=Icall+1;
            call.T{Icall}=temp.T;
            call.F{Icall}=temp.F;
            %call.score(Icall)=temp.score;
            %call.matrix=detect.matrix;
            call.ctime(Icall)=cstart+min(temp.T);
         elseif Idebug>0,
             disp(sprintf('Segment length: %6.2f, required min length: %6.2f s max length: %6.2f, FAILURE',max(temp.T)-min(temp.T), connect.seg_total_duration, ...
                 connect.seg_max_duration));
        end
       
    end %I
    %disp(sprintf('%i segments left',Icall));
elseif ~isempty(detect.T)  %Only one segment, so pass it through if it is long enough
    if max(detect.T{1})-min(detect.T{1})>=connect.seg_total_duration&max(detect.T{1})-min(detect.T{1})<=connect.seg_max_duration;
        call.F=detect.F;
        call.T=detect.T;
       % call.score=detect.score;
        %call.matrix=detect.matrix;
        call.ctime=cstart+min(detect.T{1});
        return
    end

end
if ~exist('call'),
call=[];
end



if Idebug==1,keyboard;end