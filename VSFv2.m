function [beta,varargout] = VSFv2(m,x,varargin)

% The function VSFv2 (version 2) compute the angular scattering efficiency, total 
% attenuation efficiency, total scattering efficiency and total backscattering efficiency
% for spherical or coated spherical particles.
% The angular definition could be user-defined or use pre-defined.
% The predefined angles are (in Matlab language):
% ang = [0:0.1:1 1:1:10 10:2.5:170 170:1:180] for size parameter up to 6500.
% You should define your own angles if you want a different angular resolution or the size
% parameter of the particle is greater than 6500.
%
% Usage:
%	1) beta= or [beta Qb]= or [beta Qb Qbb]= or [beta Qb Qbb Qe]=  VSFv2(m,x) 
%   2) beta= or [beta Qb]= or [beta Qb Qbb]= or [beta Qb Qbb Qe]=  VSFv2(m,x,ang)
%   3) beta= or [beta Qb]= or [beta Qb Qbb]= or [beta Qb Qbb Qe]=  VSFv2(m,x,m2,y)
%   4) beta= or [beta Qb]= or [beta Qb Qbb]= or [beta Qb Qbb Qe]=  VSFv2(m,x,m2,y,ang)
%
% Where:
%   beta = angular scattering efficiency,
%   Qb = scattering efficiency
%   Qbb = backscattering efficiency
%   Qe = attenuation efficiency
%   m = relative refractive index of the (core) particle
%   x = size parameter of the (core) particle, could be a vector
%   m2 = relative refractive index of the coating
%   y = size parameter of the ENTIRE particle (core size + coating thickness), could be a vector
%   ang = angles (in degree) at which you want your output of beta. Otherwise use predefined angular range.
%
% Note:
%   1) x could be a vector, and the beta will be a matrix.
%      if x is a vector, y has to be a vector of same length
%   2) If use the default angular range (usage 1 and 3), the first row of "beta" stores the 
%      default angle definition.
%
% What's new in version 2.
%   * use an analytical formula to calculate the backscattering efficiency, the result, however,
%     is not different from the previous calculation using integration.
%
% Written by
%		Xiaodong Zhang
%		Dalhousie University
%		January 5, 1999.


if nargin == 2
   coat = 0;
   userangle = 0;
   y = x;
elseif nargin == 4
   coat = 1;
   userangle = 0;
   m2 = varargin{1};
   y = varargin{2};
elseif nargin == 3
   coat = 0;
   userangle = 1;
   angle = varargin{1};
   y = x;
elseif nargin == 5
   coat = 1;
   userangle = 1;
   m2 = varargin{1};
   y = varargin{2};
   angle = varargin{3};
else
   error('The wrong input arguments');
end

if nargout == 1
    SCA = 0;
    BSC = 0;
    ATT = 0;
elseif nargout == 2
    SCA = 1;
    BSC = 0;
    ATT = 0;
elseif nargout == 3
    SCA = 1;
    BSC = 1;
    ATT = 0;
elseif nargout == 4
    SCA = 1;
    BSC = 1;
    ATT = 1;
else
   error('Wrong number in output arguments');
end

x = x(:);
y = y(:);

if length(x) ~= length(y)
   error('For coated case, the sizes of radius vectors (core and core+coat) should be same');
end


% Al_6500_0_180 is the Legendre coefficients calculated for size = 6500 and
% angle from 0 to 180.
% It includes:	ang0 = [0:0.1:1 1:1:10 10:2.5:170 170:1:180];
%               ang0r = ang0*pi/180;
%				lim=6500
%				nc0=Nstop(lim)
%				p0 and t0, [p0,t0]=Alegendre(ang,nc0)

if ~userangle
   load A_6500_0_180;
   if max(y) > lim
      ang0 = [0:0.1:1 1:1:10 10:2.5:170 170:1:180];
      ang0 = unique(ang0);
      ang0r = ang0*pi/180;
      nc0 = Nstop(max(y));
      [p0,t0] = ALegendr(ang0r,nc0);
   else
      lim = max(y);
      nc0 = Nstop(lim);
      p0=p0(1:nc0,:);
      t0=t0(1:nc0,:);
   end
else
   ang0 = angle(:);
   ang0 = unique(ang0');
   ang0r = ang0*pi/180;
   nc0 = Nstop(max(y));
   [p0,t0] = ALegendr(ang0r,nc0);
end

load Bb_6500
if max(y) > lim
   nc0 = Nstop(max(y));
   FE=NewBbBak(nc0);
end

beta = zeros(length(y),length(ang0));
if SCA
   Qb = zeros(size(y));
end
if BSC
   Qbb = zeros(size(y));
end
if ATT
    Qe = zeros(size(y));
end


E0 = [1:nc0]';
E10 = 2*E0+1;
E30 = E0.*(E0+1);
E20 = E10./E30; %for n=1:nc,E(n,1)=(2*n+1)/(n*(n+1)),end;

for i=1:length(y)
   if rem(i,100) == 0
      disp(['------This is No. ' num2str(i) ' calculation for radius']);
   end
   
   nc = Nstop(y(i));
   y2 = y(i)*y(i);
   if ~coat
      [a,b]=ScatCoef(m,y(i),nc);
   else
      [a,b]=ScatCoef_coat(m,m2,x(i),y(i),nc);
   end
   if SCA
      Qb(i) = E10(1:nc)'*(abs(a).^2 + abs(b).^2)*2/y2;
   end
   if ATT
       Qe(i) = E10(1:nc)'*(real(a) + real(b))*2/y2;
   end
   if BSC
      Qbb(i) = NewBbl(nc,FE,a,b)*2/y2;
%      Qbb(i) = NewBb(nc,a,b)*2/y(i)/y(i);
   end
   
   E2 = E20(1:nc);
%   p = p0(1:nc,:);
%   t = t0(1:nc,:);
   S1 = (a.*E2)'*p0(1:nc,:) + (b.*E2)'*t0(1:nc,:);
   S2 = (a.*E2)'*t0(1:nc,:) + (b.*E2)'*p0(1:nc,:);
   beta(i,:) = ((S2.*conj(S2))+(S1.*conj(S1)))/2/pi/y2;
end
if ~userangle
    beta=[ang0; beta];
end
if nargout == 2
   varargout{1} = Qb;
end
if nargout == 3
   varargout{1} = Qb;
   varargout{2} = Qbb;
end
if nargout == 4
    varargout{1} = Qb;
    varargout{2} = Qbb;
    varargout{3} = Qe;
end
return