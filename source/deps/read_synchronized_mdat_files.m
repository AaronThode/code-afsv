function [x,fparms]=read_synchronized_mdat_files(fname_in,tdate_start,tlen)
%function [x,fparms]=read_synchronized_mdat_files(fname_in,tdate_start,tlen)
% Input:
%   fname_in:  full filename and pathname
	
		
        [dirname,fname,extt] = fileparts(fname_in);
        fname=[fname extt];
        dirname=dirname(~isspace(dirname));

		thisdir=pwd;
		cd(dirname)
		[xbot,fparms_bot]=read_mdat_file(fname,tdate_start,tlen,1);
        
        fparms=fparms_bot;
        
        if isempty(fparms.synch.file)
            x=xbot;
            cd(thisdir);
            return
        end
        
		if tdate_start<0
			tdate_start=fparms_bot.synch.time;
		end
		
		elapsed_time=tdate_start-fparms_bot.synch.time;
		if elapsed_time>0
			etime=datevec(elapsed_time);
			nsec=3600*etime(4)+60*etime(5)+etime(6);
		else
			etime=datevec(-elapsed_time);
			nsec=3600*etime(4)+60*etime(5)+etime(6);
			nsec=-nsec;
		end
		
		toffset=fparms_bot.synch.offset+fparms_bot.synch.drift*nsec/(3600*1000);
		
        
		[xtop,fparms_top]=read_mdat_file(fparms.synch.file,tdate_start+datenum(0,0,0,0,0,toffset),tlen,1);
		
		x=[xbot xtop];  %Start with bottom element, move to top.
		
        rd=[fparms_bot.geom.rd fparms_top.geom.rd];
        %spacing=[fparms_bot.geom.spacing fparms_top.geom.spacing];
        Igood=[fparms_bot.Igood fparms_top.Igood];
        
        %Include only valid channels..
        rd=rd(Igood>0);
        Ichan=1:size(x,2);
        x=x(:,Igood>0);
        Ichan=Ichan(Igood>0);
        
        [rd,Iorder]=sort(rd);  %Shallowest element is first...
        %rd=fliplr(rd);
        %Igood=Igood(Iorder);
		x=x(:,Iorder);
        %Igood=find(Igood(Iorder)>0);
        spacing=[0 diff(rd)];
		
        fparms=fparms_bot;
        fparms.geom.rd=rd;
        fparms.geom.Igood=Igood;
        fparms.geom.Iorder=Iorder;
        fparms.geom.spacing=spacing;
     
		%x=flipud(x);  %Shallowest element now first element, matching head.rd
		cd(thisdir);
		
		%%Short-circuit data load with synthesized file...
		
	end
