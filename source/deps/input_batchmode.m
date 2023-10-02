function [Batch_mode, Mode_options] = input_batchmode(app, Batch_type)
%Mode_options	=	{'Process Visible Window', 'Start Bulk Processing File', 'Start Bulk Processing Folder','Load Bulk Processing'};
Mode_options	=	{'Process Visible Window', 'Start Bulk Processing', 'Load Bulk Processing'};
qstring	=	'Which operation do you want?';
title	=	Batch_type;

button	=	questdlg(qstring, title ,...
    Mode_options{1}, Mode_options{2}, Mode_options{3},...
    Mode_options{1});

Batch_mode	=	button;

switch Batch_mode
    case 'Start Bulk Processing'
        Mode_options	=	{'File','Folder'};
        title	=	Batch_type;

        button	=	questdlg(qstring, title ,...
            Mode_options{1}, Mode_options{2}, ...
            Mode_options{1});

        Batch_mode	=	[Batch_mode ' ' button];

end
end
