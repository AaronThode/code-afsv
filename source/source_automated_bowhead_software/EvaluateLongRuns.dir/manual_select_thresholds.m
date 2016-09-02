function optvec=manual_select_thresholds(feature_name);
%%%Restrict ranges and recompute successes...
for JJ=1:length(feature_name),
    prompt{(JJ-1)*2+1}=sprintf('%s min value',feature_name{JJ});
    prompt{(JJ-1)*2+2}=sprintf('%s max value',feature_name{JJ});
    defaultanswer{(JJ-1)*2+1}='0';
    defaultanswer{(JJ-1)*2+2}='1';
end
name='Feature pruning';
numlines=1;

answer=inputdlg(prompt,name,numlines,defaultanswer);

for I=1:length(answer);
    optvec(I)=str2num(answer{I});
end

end