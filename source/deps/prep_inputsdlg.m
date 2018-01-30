
function [Prompt,Formats,DefAns,Options]=prep_inputsdlg(names, Description, defaults)  %Annotation edit window setup

Options.Resize='on';
Options.WindowStyle='normal';
Options.Interpreter='tex';

Prompt = {};
Formats = {};
DefAns = struct([]);
II=1;
Iother=0;


while II<=length(names)
  %while II<=12
  
    
    Prompt(II,:) = {Description{II}, names{II},[]};
    DefAns(1).(names{II}) = defaults{II};
    
    
    switch(names{II})
        case 'author'
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'text';
            Formats(1,1).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
            %DefAns.(names{II}) = defaults{II};
%         case 'bearing'
%             Formats(1,2).type = 'edit';
%             Formats(1,2).format = 'text';
%             Formats(1,2).span = [1 1];
%             %Formats(2,1).size = 200; % automatically assign the height
%             %DefAns.(names{II}) = defaults{II};
%             
        case 'sig_type'
            Formats(1,2).type = 'list';
            Formats(1,2).format = 'text';
            Formats(1,2).style = 'radiobutton';
            Formats(1,2).items = {'FM' 'Pulsive'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.(names{II}) = 'FM';
            

        case 'call_type'
            
            Formats(2,2).type = 'edit';
            Formats(2,2).format = 'text';
            Formats(2,2).span = [1 1];
           
%             Formats(2,2).type = 'list';
%             Formats(2,2).format = 'text';
%             Formats(2,2).style = 'radiobutton';
             %Formats(2,2).items = {'upsweep' 'downsweep' 'U-shaped' 'N-shaped','other'};
             
%             % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
%             %DefAns.(names{II}) = 'downsweep';
             if ~any(strcmp(DefAns.(names{II}),Formats(2,2).items))
                 DefAns.(names{II}) = 'downsweep';
             end
            
        case 'min_freq'
            Formats(3,1).type = 'edit';
            Formats(3,1).format = 'text';
            Formats(3,1).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
            DefAns.(names{II}) = defaults{II};
            
        case 'max_freq'
            Formats(3,2).type = 'edit';
            Formats(3,2).format = 'text';
            Formats(3,2).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
            DefAns.(names{II}) = defaults{II};
            
        case 'duration'
            Formats(4,1).type = 'edit';
            Formats(4,1).format = 'text';
            Formats(4,1).span = [1 1];
            
        case 'slope'
            Formats(4,2).type = 'edit';
            Formats(4,2).format = 'text';
            Formats(4,2).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
        case 'SNR_rms_dB'
            Formats(5,1).type = 'edit';
            Formats(5,1).format = 'text';
            
            DefAns.(names{II}) =defaults{II};
        case 'num_pulses'
            Formats(6,1).type = 'list';
            Formats(6,1).format = 'text';
            %Formats(5,1).span=[1 1];
            Formats(6,1).style = 'radiobutton';
            Formats(6,1).items = {'0' '1' '2' '3'};
            DefAns.(names{II}) = '0';
            
        case 'num_harmonics'
            Formats(6,2).type = 'list';
            Formats(6,2).format = 'text';
            %Formats(5,1).span=[1 1];
            Formats(6,2).style = 'radiobutton';
            Formats(6,2).items = {'0' '1' '2' '3'};
            DefAns.(names{II}) = '0';
            
        case 'multipath'
            Formats(7,1).type = 'list';
            Formats(7,1).format = 'text';
            Formats(7,1).style = 'radiobutton';
            Formats(7,1).items = {'0' '1' '2' '3'};
            DefAns.(names{II}) = '0';
            
        case 'confidence'
            Formats(7,2).type = 'list';
            Formats(7,2).format = 'text';
            Formats(7,2).style = 'radiobutton';
            Formats(7,2).items = {'0' '1' '2' '3','4','5'};
            DefAns.(names{II}) = '5';
        
         
        case 'link_names'
            Formats(8,1).type = 'edit';
            Formats(8,1).format = 'text';
            Formats(8,1).limits = [0 3]; % multi-select files
            Formats(8,1).size = [-1 80];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(8,1).span = [1 2];  % item is 1 field x 3 fields
            
        case 'link_hashtags'
            Formats(9,2).type = 'edit';
            Formats(9,2).format = 'text';
            Formats(9,2).limits = [0 3]; % multi-select files
            Formats(9,2).size = [-1 80];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(9,2).span = [1 1];  % item is 1 field x 3 fields
            
        case 'comments'  %%May  be the last field
            
            Irow=10;
            Icol=1;
            Formats(Irow,1).type = 'edit';
            Formats(Irow,1).format = 'text';
            Formats(Irow,1).limits = [0 3]; % multi-select files
            Formats(Irow,1).size = [-1 -1];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(Irow,1).span = [1 2];  % item is 1 field x 3 fields
          
   
%         
        otherwise
            Iother=Iother+1;
%             Irow=ceil(II/2);
            Irow=ceil(Iother/2)+9;
            Icol=2-rem(Iother,2);
            Prompt(II,:) = {Description{II}, names{II},[]};
            Formats(Irow,Icol).type = 'edit';
            Formats(Irow,Icol).format = 'text';
            %Formats(Irow,Icol).limits = [0 3]; % multi-select files
            Formats(Irow,Icol).size = [-1 -1];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(Irow,Icol).span = [1 1];  % item is 1 field x 3 fields
            
    end
    II=II+1;
end

end

