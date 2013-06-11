function	Data	=	scratch(Data)

if	isempty(Data)
	return;
end

req_fields	=	{'Description', 'Template', 'Events'};
if	~isstruct(Data) || ~all(isfield(Data, req_fields))
	errmsg	=	'Existing Data must be a structure containing the fields: \n ';
	for		ii	=	1:length(req_fields)
		errmsg	=	[errmsg req_fields{ii} ', ']; %#ok<AGROW>
	end
	errmsg(end)		=	[];		errmsg(end)		=	[];
	error(errmsg, []);
end

if	~iscell(Data.Description)
	error('Description must be a cell array of strings');
end

if	~isstruct(Data.Template)
	error('Template must be a struct with the desired data field names');
end

N		=	length(Data.Description);
names	=	fieldnames(Data.Template);
if	N ~= length(names)
	error('# of fields in Template do not match # of descriptions');
end

if	~strcmp(names{1}, 'start_time')
	error('First field in Template (and Events) must be start_time');
end

if	isempty(Data.Events)
	return;
else
	names2	=	fieldnames(Data.Events);
end

if	~strcmp(names, names2)
	error('Template and Events have differring structures');
end

%	At this point assume most everything is ok
%	return

%	Code to convert all fields (except start_time to strings)
for ii	=	2:N
	%	!?!? Existing strings are automatically passed through, but do I
	%	need to specifically check for data types other than numeric?
	Data.Template.(names{ii})		=	num2str(Data.Template.(names{ii}));
	for	jj	=	1:length(Data.Events)
		Data.Events(jj).(names{ii})		=	num2str(Data.Events(jj).(names{ii}));
	end
end
