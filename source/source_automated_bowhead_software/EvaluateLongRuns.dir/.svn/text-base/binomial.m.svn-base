function p=binomial(k,n,p0)
%%probability of getting between 1 and k successes in n trials

pp=0;

%for I=1:k
%    disp(I)
p=factorial(n)./(factorial(k).*factorial(n-k));
p=p.*(p0^k).*((1-p0).^(n-k))
%    pp=pp+p;
%end

if 1==0
    pf=0;
    for k=2:n
        pf=pf+binomial(k,n,p0);
    end
end