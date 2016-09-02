function run_energy_detector(fname,param)
%write Cshell to run energy detector on fname.
% Input:
%   fname: complete file name (including path).
%   param: structure with parameter  fields immeadiately below (no
%   param.energy)  

%fprintf(f_det,'foreach i (%s/*%s)       \n',params.file_dir, params.exten);

if ~isempty(findstr(lower(computer),'mac')),
    slashstr='/';
else
    slashstr='\';
end

Islash=max(findstr(fname,slashstr));
param.exten=fname((Islash+1):end);
param.file_dir=fname(1:(Islash-1));


scriptname=writeEnergyDetectorCshell(param);
%%Run script if energy parameters change...
disp('Running JAVA script');
%!rm *snips *detsum
eval(sprintf('! ./%s > outtt.txt',scriptname));
%! ./masterEDET.scr > outtt.txt



end

