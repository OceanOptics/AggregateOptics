% logarithmic derivative
% Dn(x) = d(ln(phi(x))/dx;
% Dn-1(x) = n/x - 1/(Dn(x)+n/x)
% Written by and copyright
%     Xiaodong ZHang  zxd@predator.ocean.dal.ca
%     Dept. of Oceanography
%     Dalhousie University
%     June 1, 1997
% Modified on May 14, 1998
% nst value has been increase by 45 from original calculation.
% This increasement will increase the accuracy of dn estimate for absorbing case
% and with large size.

function dn = dlnphi(rho, nmax)
%nst = ceil(max(nmax,abs(rho)))+15; this is original equation
nst = ceil(max(nmax,abs(rho)))+50;
dn = zeros(nst,1);
dn(nst) = 0;
for n=nst:-1:2
	dn(n-1) = n/rho - 1/(dn(n)+n/rho);
end
dn0 = 1/rho - 1/(dn(1)+1/rho);
dn00 = cos(rho)/sin(rho);
fac = dn00/dn0;
dn = dn(1:nmax)*fac;