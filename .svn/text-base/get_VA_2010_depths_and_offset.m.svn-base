function [rd,D,Iorder,torigin,toffset,tdrift,Iwant]=get_VA_2010_depths_and_offset(year,myfile,myscenario)
% Iwant: channels desired, used to remove bad channels from rd and x..
%Note that Iorder(1) is deepest element, but rd(1) is shallowest element...

if ~exist('myfile','var')
    myfile=[];
end
if ~exist('myscenario','var')
    myscenario=[];
end
torigin=[];
toffset=[];
tdrift=[];

switch year
    case 2010
        %%spacing is distance between elements, starting from ocean bottom.
        spacing=.3048*([13+5/12 9+10/12 10  9+11.5/12 9+9/12 9+6/12 5+3/12 ...
            5 4+6/12 4+5/12 4+5/12 10+3/12 10+5/12 10 9+7.5/12]);  %Convert from inches to meters
        
        Iorder=[1 2 3 4 5 6 9 7 10 8 11 12 13 14 15 16];  %Element 1 is ocean bottom, 16 nearest the surface...
        Iwant=1:15;  %Iwant(1) is ocean bottom, etc.
        
        %%Original SAGA inversion before correcting *.dat file 0 m depth
        %         D=54.87;
        %         rd=fliplr(D-cumsum(spacing));
        %         rd=rd-(15.6-13.5)+0.2693;
        
        %Short-range whale MFP inversion (1.5 km)
        %         D=56.44;
        %         rd=fliplr(D-cumsum(spacing));
        %         rd=rd+(15.485-17.9336);
        
        
        %%Original results when depth was allowed to be inverted along with large changes in element depth...
        
        %         if ~isempty(findstr(lower(myscenario),'close'))||~isempty(findstr(lower(myscenario),'short'))
        %             disp('Selecting short range inversion');
        %             D=56.44;
        %             rd=fliplr(D-cumsum(spacing));
        %             rd=rd+(15.485-17.9336);
        %         elseif ~isempty(findstr(lower(myscenario),'mid'))
        %             disp('Selecting mid range inversion');
        %
        %             D=55;
        %             rd=fliplr(D-cumsum(spacing));
        %             rd=rd-(16.4936-12.143);
        %         elseif ~isempty(findstr(lower(myscenario),'far'))
        %             disp('Selecting far range inversion');
        %
        %             D=55;
        %             rd=fliplr(D-cumsum(spacing));
        %             rd=rd+(15.485-17.9336)+0.2;
        %
        %         elseif ~isempty(findstr(lower(myscenario),'ultra'))
        %             disp('Selecting ultra range inversion');
        %
        %             D=54.5;
        %             rd=fliplr(D-cumsum(spacing));
        %             rd=rd-15.9936+16.723;
        %         else
        %             disp('Selecting default  inversion');
        %
        %             D=55;
        %             rd=fliplr(D-cumsum(spacing));
        %             rd=rd-(16.4936-12.143);
        %         end
        
        D=55;
        rd=fliplr(D-cumsum(spacing));
              
        if ~isempty(findstr(lower(myscenario),'close'))||~isempty(findstr(lower(myscenario),'short'))
            disp('Selecting short range inversion');
            rd=rd-min(rd)+14.143;
        elseif ~isempty(findstr(lower(myscenario),'mid'))
            disp('Selecting mid range inversion');
            rd=rd-min(rd)+13;
        elseif ~isempty(findstr(lower(myscenario),'far'))
            disp('Selecting far range inversion');
             rd=rd-min(rd)+13.857;
        elseif ~isempty(findstr(lower(myscenario),'ultra'))
            disp('Selecting ultra range inversion');
            rd=rd-min(rd)+14.714;
%         else
%             disp('Selecting default  inversion');
%             rd=rd-(16.4936-12.143);
        end
        
       
end

if findstr(myfile,'103.MDAT')
    torigin=datenum(2010,8,31,10,00,29);
    toffset=3-0.31648-0.17888+0.00132; %20100831T100029
    %toffset=3-0.31648-0.17888+0.01744-0.00016; %soundsamp_20100831T114524
    
    tdrift=17.3*3600/6295; %msec per hour
    
elseif findstr(myfile,'44.MDAT')
    torigin=datenum(2010,8,20,1,42,38);
    toffset=0.5-0.15-0.005-0.5808/1500+2*0.0001655-1/6250;
    tdrift=0;
    
    Iwant(3)=[];
    rd(end-2)=[];
elseif findstr(myfile,'104.MDAT')
    
    torigin=datenum(2010,8,31,14,58,54);
    toffset=2+0.5-0.13-0.011-1/6250;
    tdrift=0;
    
    Iwant(3)=[];
    rd(end-2)=[];
elseif findstr(myfile,'101.MDAT')
    
    torigin=datenum('31-Aug-2010 00:47:50');
    toffset=2+.2269+0.0091+7/6250;
    tdrift=0;
    
    rd(13)=[];
    rd([11 13])=[];
    Iwant([2 3 5])=[];
    
    
    
    
end