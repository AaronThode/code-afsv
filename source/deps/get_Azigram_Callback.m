function Batch_vars=get_Azigram_Callback(handles)

Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
Batch_desc{2}='Bearing Range Color Scale';
Batch_desc{3}='Bearing bias/correction (default is 11.7 degrees)';
Batch_desc{4}='Bearing algorithm multiplicative (Default 1)';

Batch_vars.sec_avg	=	'0.1';
Batch_vars.climm='[0 360]';
Batch_vars.brefa='11.7';
handBatch_vars.alg='1';

if isfield(handles,'azigram')
    if isfield(handles.azigram,'sec_avg')
        Batch_vars.sec_avg=handles.azigram.sec_avg;
    end
    if isfield(handles.azigram,'climm')
        Batch_vars.climm=handles.azigram.climm;
    end
    
    if isfield(handles.azigram,'brefa')
        Batch_vars.brefa=handles.azigram.brefa;
    end
    
     if isfield(handles.azigram,'alg')
        Batch_vars.alg=handles.azigram.alg;
    end
end

Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, 'Vector Sensor Processing');
end