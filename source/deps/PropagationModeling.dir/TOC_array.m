%%%%%%%%%%%TOC_array.m%%%%%%%%%%
%  Select a source array--currently hardwired to 'single_element'
%function [x,y,Sf_pos,X,Y]=TOC_array(case_array),

function [x,y,Sf_pos,X,Y]=TOC_array(case_array),
      
if strcmp(case_array,'select_from_menu'),
    menu_chc={'single_element','Fairfield_1680','Ewing_20'};
    Ichc=menu('Select an airgun array profile:',menu_chc)
    case_array=menu_chc{Ichc};
end

switch case_array,
    case{'single_element'},
        y=[0];
        x=[0 ];
        Sf_pos=1;
        x=x-mean(x);
        y=y-mean(y);
        [X,Y]=meshgrid(x,y);
    case{'Fairfield_1680'},
        y=[-2.5 2.5];
        x=[0 3 6 9];
        Volume=[160 360 80 240;
            160 360 80 240]/80;
        Sf_pos=Volume.^(1/3);  % Individual source strengths, relative to smallest element, 
        x=x-mean(x);
        y=y-mean(y);
        [X,Y]=meshgrid(x,y);
        % scaled by cube of volume
    case{'Ewing_20'},
        %20 gun ewing deployment...
        %         raw=[145 -16.8 -4.6;
        %             850 -15.2 4.6;
        %             385 -13.7 0;
        %             235 -12.2 -4.6;
        %             640 -10.7  4.6;
        %             585 -9.1   0;
        %             250 -7.6  -4.6;
        %             850 -6.1   4.6;
        %             540 -3.9   0.0;
        %             250 -2.4 -4.6;
        %             235  2.4 -4.6;
        %             585  3.9  0.0;
        %             850  6.1  4.6;
        %             120  7.6 -4.6;
        %             500  9.1  0.0;
        %             80   10.7 4.6;
        %             200  12.2 -4.6;
        %             305  13.7  0.0;
        %             850  15.2  4.6;
        %             145  16.8 -4.6;];
        raw=[ 210785000000.00 145 -16.8 -4.6
            368360000000.00 850 -15.2 4.6
            326340000000.00 385 -13.7 0
            263310000000.00 235 -12.2 -4.6
            363875652173.91 640 -10.7 4.6
            359355714285.71 585 -9.1 0
            270472500000.00 250 -7.6 -4.6
            368360000000.00 850 -6.1 4.6
            354444285714.29 540 -3.9 0
            270472500000.00 250 -2.4 -4.6
            263310000000.00 235 2.4 -4.6
            359355714285.71 585 3.9 0
            368360000000.00 850 6.1 4.6
            189365714285.71 120 7.6 -4.6
            349896666666.67 500 9.1 0
            149856000000.00 80 10.7 4.6
            245738000000.00 200 12.2 -4.6
            294688571428.57 305 13.7 0
            368360000000.00 850 15.2 4.6
            210785000000.00 145 16.8 -4.6];
        
        Volume=raw(:,2)/80;
        Y=raw(:,3);
        X=raw(:,4);
        I80=find(Volume==1);
        Sf_pos_0=Volume.^(1/3);  % Individual source strengths, relative to smallest element, 
        Sf_pos=raw(:,1)/raw(I80,1);  %% Strength relative to the 80 cu in element
        % scaled by cube of volume
    case{'Kondor_3930'},
        
        raw=[300 0   0;
            300  3   0;
            100  5   0;
            40   7   0;
            70   9   0;
            200 11   0;
            500  14   0;
            
            500 0   12.5;
            90 3  12.5;
            60 5  12.5;
            150 7 12.5;
            40  9  12.5;
            70  11  12.5;
            250 14 12.5;
            
            
            300 0   25;
            300  3   25;
            70  5   25;
            70   7   25;
            70   9   25;
            200 11   25;
            500  14   25;];
        Volume=raw(:,1)/80;
        X=raw(:,2);
        Y=raw(:,3);
        Sf_pos=Volume.^(1/3);  % Individual source strengths, relative to smallest element, 
        % scaled by cube of volume
end