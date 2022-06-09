% RB1	the Ricatti-Bessel function of the first kind
%		RB1(rho, nmax) for the value rho from n=1 to n=nmax.
% Written by and copyright
%		Dave Barnett
%		Optical Engineering Group
%		Dept. of Electrical and Electronic Engineering
%		Loughborough University
%		20th November 1996

function phi = RB(rho, nmax)
phi = zeros(nmax,1);
phi_1 = cos(rho);
phi_0 = sin(rho);
phi(1) = phi_0/rho - phi_1;
phi(2) = 3*phi(1)/rho - phi_0;
for n=2:(nmax-1)
	phi(n+1) = (2*n+1)*phi(n)/rho - phi(n-1);
end


