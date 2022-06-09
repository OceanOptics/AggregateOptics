function [a,b] = ScatCoef_coat(m1,m2,x,y,nmax)
% ScatCoef_coat	Scattering coefficients for coated sphere
%		[a,b]=ScatCoef_coat(m1,m2,x,y,nmax)
%		returns the two column vectors a and b containing
%		the scattering coefficients for coated particles of 
%     core size x and coat size y as well as core refractive
%     indes m1 and coat index m2, from n=1 to n=nmax.
% The biggest improvement:
%		1) calculate chi using stable upward recursion;
%		2) calculate Dn using stable downward recursion;
%		3) calcualte phi using 1/phi = chi*Dn - dchi,
%		4) where Dn is derivative of ln(phi) and dchi derivative of chi.
% Written by and copyright
%		Xiaodong Zhang
%		Dept. of Oceanography
%     Dalhousie University
%     May 13, 1998
% Ref. Absorption and Scattering of Light by Small Particles
%     by Craig F. Bohren and Donald R. Huffman
%     [a,b]=ScatCoef_coat(m1,m2,x,y,nmax,method)
%     


if nargin  ~= 5
   error('Require 5 input arguments.')
end

m = m2/m1;

chiy = RB2(y, nmax);
chim2y = RB2(m2*y, nmax);
chim2x = RB2(m2*x, nmax);

dchiy = DRB(chiy,y,2);
dchim2y = DRB(chim2y,m2*y,2);
dchim2x = DRB(chim2x,m2*x,2);

%dny = dlnphi_u(y,nmax);
%dnm2y = dlnphi_u(m2*y,nmax);
%dnm1x = dlnphi_u(m1*x,nmax);
%dnm2x = dlnphi_u(m2*x,nmax);

dny = dlnphi(y,nmax);
dnm2y = dlnphi(m2*y,nmax);
dnm1x = dlnphi(m1*x,nmax);
dnm2x = dlnphi(m2*x,nmax);

%phiy = RB(y,nmax);
phim2x = RB(m2*x,nmax);
phim2y = RB(m2*y,nmax);

tmpa = chim2x.*dnm2x - dchim2x;
tmpb = chim2y.*dnm2y - dchim2y;
inda = find(tmpa ~= 0);
indb = find(tmpb ~= 0);

phiy = 1./(chiy.*dny - dchiy);
phim2x(inda) = 1./tmpa(inda);
phim2y(indb) = 1./tmpb(indb);

xiy = phiy - i*chiy;

A = phim2x.*(m*dnm1x-dnm2x) ./ (m*dnm1x.*chim2x-dchim2x);
B = phim2x.*(m*dnm2x-dnm1x) ./ (m*dchim2x-dnm1x.*chim2x);
   
D = (dnm2y-A.*dchim2y./phim2y) ./ (1-A.*chim2y./phim2y);
G = (dnm2y-B.*dchim2y./phim2y) ./ (1-B.*chim2y./phim2y);
   
tmp = [1:nmax]';
tmpa = D/m2 + tmp/y;  
tmpb = m2*G + tmp/y;
phiyn_1 = shift_r(phiy,y,1);
xiyn_1 = shift_r(xiy,y,3);
   
a = (tmpa.*phiy - phiyn_1) ./ (tmpa.*xiy - xiyn_1);
b = (tmpb.*phiy - phiyn_1) ./ (tmpb.*xiy - xiyn_1);
return
