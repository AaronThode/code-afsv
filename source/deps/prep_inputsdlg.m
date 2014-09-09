Options.Resize='on';
Options.WindowStyle='normal';
Options.Interpreter='tex';

Prompt = {};
Formats = {};
DefAns = struct([]);
II=1;

Prompt(II,:) = {Description{II}, 'author',[]};
Formats(1,1).type = 'edit';
Formats(1,1).format = 'text';
Formats(1,1).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns(1).(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II},'sig_type',[]};
Formats(2,1).type = 'list';
Formats(2,1).format = 'text';
Formats(2,1).style = 'radiobutton';
Formats(2,1).items = {'FM' 'Pulsive'};
% Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'call_type',[]};
Formats(2,2).type = 'list';
Formats(2,2).format = 'text';
Formats(2,2).style = 'radiobutton';
Formats(2,2).items = {'upsweep' 'downsweep' 'U-shaped' 'N-shaped'};
% Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
DefAns.(names{II}) = 'downsweep';

II=II+1;
Prompt(II,:) = {Description{II}, 'min_freq',[]};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'text';
Formats(3,1).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'max_freq',[]};
Formats(3,2).type = 'edit';
Formats(3,2).format = 'text';
Formats(3,2).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'duration',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'text';
Formats(4,1).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'slope',[]};
Formats(4,2).type = 'edit';
Formats(4,2).format = 'text';
Formats(4,2).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'num_pulses',[]};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'text';
Formats(5,1).span = [1 1];
%Formats(2,1).size = 200; % automatically assign the height
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'num_harmonics',[]};
Formats(5,2).type = 'list';
Formats(5,2).format = 'text';
Formats(5,2).style = 'radiobutton';
Formats(5,2).items = {'0' '1' '2' '3','4','5','6'};
DefAns.(names{II}) = '0';

II=II+1;
Prompt(II,:) = {Description{II}, 'multipath',[]};
Formats(6,1).type = 'list';
Formats(6,1).format = 'text';
Formats(6,1).style = 'radiobutton';
Formats(6,1).items = {'0' '1' '2' '3'};
DefAns.(names{II}) = '0';

II=II+1;
Prompt(II,:) = {Description{II}, 'confidence',[]};
Formats(6,2).type = 'edit';
Formats(6,2).format = 'text';
Formats(6,2).style = 'radiobutton';
Formats(6,2).items = {'0' '1' '2' '3','4','5'};
DefAns.(names{II}) = '5';


II=II+1;
Prompt(II,:) = {Description{II}, 'link_names',[]};
Formats(7,1).type = 'edit';
Formats(7,1).format = 'text';
Formats(7,1).limits = [0 3]; % multi-select files
Formats(7,1).size = [-1 80];
%Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
Formats(7,1).span = [1 2];  % item is 1 field x 3 fields
DefAns.(names{II}) = defaults{II};

II=II+1;
Prompt(II,:) = {Description{II}, 'link_hashtags',[]};
Formats(8,1).type = 'edit';
Formats(8,1).format = 'text';
Formats(8,1).limits = [0 3]; % multi-select files
Formats(8,1).size = [-1 80];
%Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
Formats(8,1).span = [1 2];  % item is 1 field x 3 fields
DefAns.(names{II}) = defaults{II};

% We will assume that remaining fields are simple edit text boxes
Icount=0;
for II=(II+1):length(names)
    Icount=Icount+1;
    Irow=8+ceil(Icount/2);
    Icol=2-rem(Icount,2);
    Prompt(II,:) = {Description{II}, names{II},[]};
    Formats(Irow,Icol).type = 'edit';
    Formats(Irow,Icol).format = 'text';
    %Formats(Irow,Icol).limits = [0 3]; % multi-select files
    Formats(Irow,Icol).size = [-1 -1];
    %Formats(7,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
    Formats(9,1).span = [1 1];  % item is 1 field x 3 fields
    DefAns.(names{II}) = defaults{II};
end



