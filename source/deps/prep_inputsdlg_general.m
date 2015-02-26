function prep_inputsdlg_general(Description,names,defaults)

Options.Resize='on';
Options.WindowStyle='normal';
Options.Interpreter='tex';

Prompt = {};
Formats = {};
DefAns = struct([]);
II=1;

while II<=length(names)
    
    
    Prompt(II,:) = {Description{II}, names{II},[]};
    DefAns(1).(names{II}) = defaults{II};
    
    
    switch(names{II})
        case 'author'
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'text';
            Formats(1,1).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
            %DefAns.(names{II}) = defaults{II};
            
        case 'sig_type'
            Formats(2,1).type = 'list';
            Formats(2,1).format = 'text';
            Formats(2,1).style = 'radiobutton';
            Formats(2,1).items = {'FM' 'Pulsive'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            
        case 'call_type'
            Formats(2,2).type = 'list';
            Formats(2,2).format = 'text';
            Formats(2,2).style = 'radiobutton';
            Formats(2,2).items = {'upsweep' 'downsweep' 'U-shaped' 'N-shaped','other'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            %DefAns.(names{II}) = 'downsweep';
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
            
        case 'num_pulses'
            Formats(5,1).type = 'edit';
            Formats(5,1).format = 'text';
            Formats(5,1).span = [1 1];
            %Formats(2,1).size = 200; % automatically assign the height
            
        case 'num_harmonics'
            Formats(5,2).type = 'list';
            Formats(5,2).format = 'text';
            Formats(5,2).style = 'radiobutton';
            Formats(5,2).items = {'0' '1' '2' '3','4','5','6'};
            DefAns.(names{II}) = '0';
            
        case 'multipath'
            Formats(6,1).type = 'list';
            Formats(6,1).format = 'text';
            Formats(6,1).style = 'radiobutton';
            Formats(6,1).items = {'0' '1' '2' '3'};
            DefAns.(names{II}) = '0';
            
        case 'confidence'
            Formats(6,2).type = 'edit';
            Formats(6,2).format = 'text';
            Formats(6,2).style = 'radiobutton';
            Formats(6,2).items = {'0' '1' '2' '3','4','5'};
            DefAns.(names{II}) = '5';
            
        case 'link_names'
            Formats(7,1).type = 'edit';
            Formats(7,1).format = 'text';
            Formats(7,1).limits = [0 3]; % multi-select files
            Formats(7,1).size = [-1 80];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(7,1).span = [1 2];  % item is 1 field x 3 fields
            
        case 'link_hashtags'
            Formats(8,1).type = 'edit';
            Formats(8,1).format = 'text';
            Formats(8,1).limits = [0 3]; % multi-select files
            Formats(8,1).size = [-1 80];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(8,1).span = [1 2];  % item is 1 field x 3 fields
            
        case 'comments'
            Formats(9,1).type = 'edit';
            Formats(9,1).format = 'text';
            Formats(9,1).limits = [0 3]; % multi-select files
            Formats(9,1).size = [-1 -1];
            %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
            Formats(9,1).span = [1 1];  % item is 1 field x 3 fields
            
        otherwise
           
            Irow=ceil(II/2);
            Icol=2-rem(II,2);
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


