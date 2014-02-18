function close_all_figures
figno=get(0,'child');
figno=figno(figno<150);
for I=1:length(figno)
    close(figno(I));
end