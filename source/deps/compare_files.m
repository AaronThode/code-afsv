%%compare_files.m
% Compare files in this folder to JAVA_IO and other directories to ensure
% consistency in application...
%  run inside deps folders for AFSV

clear all
fnames=dir('*.m');


test_dir_all={'../../JAVA_IO/','../../File_IO/','../../UsefulScripts.dir/','../../ICI_detectors.dir','../../CreakDetector.dir'};

for J=1:length(test_dir_all)
    test_dir=test_dir_all{J};
    for I=1:length(fnames)
        fcopy=dir([test_dir fnames(I).name]);
        if length(fcopy)>0
            %disp(fcopy.name);
            eval_str=sprintf('!diff %s %s > diffout.txt',fnames(I).name,[test_dir fcopy.name]);
            eval(eval_str);
            head=dir('diffout.txt');
            fprintf('%s difference is %i bytes\n',fcopy.name,head.bytes);
            if head.bytes>0
                edit diffout.txt
                
                yes=menu(sprintf('Replace %s in %s?',fcopy.name,test_dir),'Yes','No');
                if yes==1
                    eval(sprintf('!cp %s %s',fnames(I).name,[test_dir]));
                end
                
            end
            !rm diffout.txt
            
            
        end
        
    end
end