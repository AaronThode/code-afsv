%%%%%%%alter_parameters.m%%%%
function [param,change_flag]=alter_parameters(param)

names=fieldnames(param);
names{length(names)+1}='Exit';
Ichc=menu('Select a parameter:',names);
change_flag=0;
double_level=0;
if ~strcmp(names{Ichc},'Exit')&&isstruct(param.(names{Ichc}))
    double_level=1;
    names2=fieldnames(param.(names{Ichc}));
    Ichc2=menu('Select a parameter:',names2);
else
    double_level=0;

end


while Ichc~=length(names),
    change_flag=1;

    if double_level==0,
        current_val=param.(names{Ichc});
    else
        current_val=param.(names{Ichc}).(names2{Ichc2});
    end

    for I=1:length(current_val),

        if double_level==0,
            val=input(sprintf('Change %s = %6.2f into what? (Return to leave unchanged)',names{Ichc},current_val(I)));
            if ~isempty(val)
                param.(names{Ichc})(I)=val;
            end

        else
            val=input(sprintf('Change %s.%s = %6.2f into what? (Return to leave unchanged)',names{Ichc},names2{Ichc2},current_val(I)));
            if ~isempty(val)
                param.(names{Ichc}).(names2{Ichc2})(I)=val;
            end


        end
    end
    Ichc=menu('Select a parameter:',names);

    double_level=0;
    if ~strcmp(names{Ichc},'Exit')&&isstruct(param.(names{Ichc}))
        double_level=1;
        names2=fieldnames(param.(names{Ichc}))
        Ichc2=menu('Select a parameter:',names2);
    else
        double_level=0;

    end
end

end