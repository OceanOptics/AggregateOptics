function [a,b] = ScatCoef(m,x,nmax)

% ScatCoef	Scattering coefficients
%		[a,b]=ScatCoef(m,x,nmax)
%		returns the two column vectors a and b containing
%		the scattering coefficients for particles of size x
%		and refractive index relative to medium m, from n=1
%		to n=nmax.
% Written by Xiaodong ZHang
%     Ref. Absorption and Scattering of Light by Small Particles
%     by Craig F. Bohren and Donald R. Huffman
%     [a,b]=ScatCoef(m,x,nmax,method)
%     June 29, 1997

if nargin ~= 3
   error('[a,b]=ScatCoef(m,x,nmax)')
end

phi = RB1(x, nmax);
chi = RB2(x, nmax);
xi = phi - i*chi;
dnm=dlnphi(m*x,nmax);

tmp = [1:nmax]';
tmpa = dnm/m + tmp/x;
tmpb = dnm*m + tmp/x;
phin_1 = shift_r(phi,x,1);
xin_1 = shift_r(xi,x,3);
a = (tmpa.*phi - phin_1) ./ (tmpa.*xi - xin_1);
b = (tmpb.*phi - phin_1) ./ (tmpb.*xi - xin_1);
return

