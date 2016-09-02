%%%%hausdorf_test.m%%%%
clear; close all
load haus_test

X=initial_final_image(:,1:140);
%X(:,1:90)=0;
Z=initial_final_image(:,141:(end-10));

Nshifts=size(X,2);

%[dist]=similarity_function(X, X);


[dist2]=similarity_function(fliplr(X), X);

[dist3]=similarity_function(Z, X);

%LowF=circshift(X,[10 0]);
%[dist4]=similarity_function(LowF, X);

Partial=Z(:,1:40);
[dist4]=similarity_function(Partial, X);

[dist4]=similarity_function(fliplr(Partial), X);
