function [HH]=hausdorff_image(Y, X,Iplot)

HH=Inf;
if sum(sum((X)))==0 | sum(sum((Y)))==0
    return
end
Mcol=max([size(Y,2) size(X,2)]);

if size(Y,2)~=Mcol,
    Y=[Y zeros(size(X,1),Mcol-size(Y,2))];
    %Y=logical(Y);
else
    X=[X zeros(size(X,1),Mcol-size(X,2))];
   % X=logical(X);
end


[dist{1},shifft{1}]=hausdorff_translate_mean(X, Y);
[dist{2},shifft{2}]=hausdorff_translate_mean(Y, X);
HH=min([min(dist{1}) min(dist{2})]);

if Iplot>0
    figure;
    subplot(3,1,1)
    imagesc(X);
    subplot(3,1,2)
    imagesc(Y)
    subplot(3,1,3)
    plot(shifft{1},dist{1},shifft{2},dist{2});grid on
    ylim([0 2])
    disp(sprintf('Minimum value of both hausdorff transforms is %6.2f',min([min(dist{1}) min(dist{2})])));
    pause;
    close
end