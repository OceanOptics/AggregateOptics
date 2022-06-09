% ALegendr	the angular dependent Associated Legendre Polynomials
%		[p,t]=ALegendr(ang, nmax)
%		produces matrices p and t with rows n=1 to n=nmax
%		for pi and tau functions rescpectively.


function [p,t] = ALegendr(ang, nmax)
p=zeros(length(1:nmax),length(ang));
t=p;
u=cos(ang);
p(1,:) = 1;
t(1,:) = u;
p(2,:) = 3*u;
t(2,:) = 2*u.*p(2,:)-3;
for n=3:nmax
	p(n,:) = ((2*n-1)*u.*p(n-1,:) - n*p(n-2,:))/(n-1);
	t(n,:) = n*u.*p(n,:) - (n+1)*p(n-1,:);
end
