function Batch_vars=get_Azigram_Callback(param)

%%%allow either manual input of parameters or extract parameters from data
%%%(e.g. brefa from an orientation file)
Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
Batch_desc{2}='Bearing Range Color Scale';
Batch_desc{3}='Bearing bias/correction ';
Batch_desc{4}='Bearing algorithm multiplicative (Default 1)';
Batch_desc{5}='Bearing Mask?';

Batch_vars.sec_avg	=	'0.1';
Batch_vars.climm='[0 360]';
Batch_vars.brefa='0';
Batch_vars.alg='1';
Batch_vars.mask=false;

if isfield(param,'sec_avg')
    Batch_vars.sec_avg=param.sec_avg;
end
if isfield(param,'climm')
    Batch_vars.climm=param.climm;
end

if isfield(param,'brefa')
    Batch_vars.brefa=param.brefa;
    if isfield(param,'myfile')
        if contains(param.myfile,'drifter')
            [brefa,prefa]=get_drifter_orientation(param.myfile);
            Batch_vars.brefa=brefa;
        %keyboard
        end
    end
else
    Batch_vars.brefa='0';
end

if isfield(param,'alg')
    Batch_vars.alg=param.alg;
end

if isfield(param,'mask')
    Batch_vars.mask=param.mask;
end


Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, 'Vector Sensor Processing');
end