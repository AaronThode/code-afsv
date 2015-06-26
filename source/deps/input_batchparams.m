% Asks user for specified batch processing parameters-return strings
function	[Batch_vars, Batch_names]	=	input_batchparams(Batch_vars, Batch_desc, Batch_type, num_out)

if	isstruct(Batch_vars)
    names		=	fieldnames(Batch_vars);
elseif	iscell(Batch_vars)
    names		=	Batch_vars;
    Batch_vars	=	struct;
else
    error('Batch_vars must be either a structure, or cell array of variable names');
end

N_fields	=	length(names);
if	length(Batch_desc) ~= N_fields
    error('# of Descriptors differs from # of fields');
end

%	Default values
defaults	=	cell(N_fields,1);
if	isstruct(Batch_vars)
    for	ii	=	1:N_fields
        value	=	Batch_vars.(names{ii});
        defaults{ii}	=	num2str(value);
    end
end

%	Size of each prompt field
num_lines		=	1;

%	Title for window
dlgTitle	=	['Parameters for  ' Batch_type];


%	Create input dialogbox
%answer		=	inputdlg(Batch_desc, dlgTitle, num_lines, defaults,'on');
try
    %prep_inputsdlg;
    Options.Resize='on';
    Options.WindowStyle='normal';
    Options.Interpreter='tex';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    Nrows_max=12;
    
    
    for II=1:length(names)
        
        Prompt(II,:) = {Batch_desc{II}, names{II},[]};
        DefAns(1).(names{II}) = defaults{II};
        
        Icol=1;
        Irow=II;
        if II>Nrows_max
            Icol=2;
            Irow=II-Nrows_max;
            
        end
        Formats(Irow,Icol).type = 'edit';
        Formats(Irow,Icol).format = 'text';
        Formats(Irow,Icol).span = [1 1];
        
    end
    [answer,Cancelled] = inputsdlg(Prompt,dlgTitle,Formats,DefAns,Options);
    
    %	Put variables into output structure
    if	isempty(answer)
        Batch_vars	=	[];
        return;
    end
    
    for	ii	=	1:N_fields
        Batch_vars.(names{ii})	=	answer.(names{ii});
        if nargin>3
            try
                Batch_vars.(names{ii})	=	eval(answer.(names{ii}));
            end
            
        end
    end
catch
    answer		=	inputdlg(Batch_desc, dlgTitle, num_lines, defaults,'on');
    for	ii	=	1:N_fields
        Batch_vars.(names{ii})	=	answer{ii};
        if nargin>3
            try
                Batch_vars.(names{ii})	=	eval(answer{ii});
            end
            
        end
    end
    
end



Batch_names=fieldnames(Batch_vars);

end
