%This loads original UNIX path
path1=getenv('PATH');
path1=sprintf('%s:%s','/usr/local/bin',path1);

path1=sprintf('%s:%s',path1,'/Users/thode/Projects/PropagationCodes/at/UnixScripts');
path1=sprintf('%s:%s',path1,'/Users/thode/Projects/PropagationCodes/at/bin');
path1=sprintf('%s:%s',path1,'/Users/thode/Projects/PropagationCodes/RAM');

setenv('PATH',path1)
%setenv('PATH', [getenv('PATH') ';D:\Perl\bin']);


%%Allow gfortran compiled programs to run within MATLAB shell
%http://www.mathworks.com/matlabcentral/answers/12822-shell-command-on-matlab-command-line-shows-error
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');

%  From: http://www.mathworks.co.uk/matlabcentral/answers/44388-or-system-or-unix-input-redirection
%setenv('GFORTRAN_STDIN_UNIT', '5') ;
%setenv('GFORTRAN_STDOUT_UNIT', '6') ;
%setenv('GFORTRAN_STDERR_UNIT', '0');

%setenv(?GFORTRAN_STDIN_UNIT?, ?-1?)
%setenv(?GFORTRAN_STDOUT_UNIT?, ?-1?)
%setenv(?GFORTRAN_STDERR_UNIT?, ?-1?)
