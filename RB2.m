% RB2	the Ricatti-Bessel function of the second kind
%		RB2(rho, nmax) for the value rho from n=1 to n=nmax.
% Written by and copyright
%		Dave Barnett
%		Optical Engineering Group
%		Dept. of Electrical and Electronic Engineering
%		Loughborough University
%		20th November 1996

function zeta = RB2(rho, nmax)
zeta = zeros(nmax,1);
zeta_1 = -sin(rho);
zeta_0 = cos(rho);
zeta(1) = zeta_0/rho - zeta_1;
zeta(2) = 3*zeta(1)/rho - zeta_0;
for n=2:(nmax-1)
	zeta(n+1) = (2*n+1)*zeta(n)/rho - zeta(n-1);
end