% calculate the differential Raccatti Bessel function using recursion
% drb(n) = rb(n-1) - n*rb(n)/x
% rb(0) = rb0 = sin(x) for phi (Bessel of first order)
%             = -cos(x) for zeta (Bessel of second order)
%             = sin(x) - cos(x)i for xi (Bessel of third order)

% Xiaodong Zhang
% Dept. of Oceanography
% Dalhousie University
% July 4, 1997

function drb=DRB(rb,x,order)

if nargin < 3
   error('There must be three input arguments');
elseif nargin == 3
   rbn_1 = shift_r(rb,x,order);
else
   error('Too much input arguments');
end

rbn_1 = shift_r(rb,x,order);
en = [1:length(rb)]';

drb = rbn_1 - en.*rb/x;


