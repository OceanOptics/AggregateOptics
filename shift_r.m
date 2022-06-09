function rbn_1=shift_r(rb,x,order)
% Make the Bessel series of rb one position prior.
% Xiaodong Zhang
% The order is 1, 2 or 3 defining the order of Bessel function.

rbn_1=rb;
len = length(rb);

if order == 1
   rbn_1(1) = sin(x);
elseif order == 2
   rbn_1(1) = cos(x);
elseif order == 3
   rbn_1(1) = sin(x)-i*cos(x);
else
   error('The wrong input for "order"')
end

rbn_1(2:len) = rb(1:len-1);
