function data_out=merge_feature_data(data_out,data)

names=fieldnames(data);
for I=1:length(names)
    if ~isfield(data_out,names{I})
        data_out.(names{I})=[data.(names{I})];
    else

        data_out.(names{I})=[data_out.(names{I}) data.(names{I})];
    end
end