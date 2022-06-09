function FE=NewBbBak(nc0)
% Pre-calculate the coefficients needed for the NewBb and save them as
% a matfile for later use.
% Xiaodong Zhang
% May 5, 1999

E = [1:nc0]';

FE = dfractorial(E);

return

function y=dfractorial(x)
y = x;
for i=1:length(x)
   if rem(x(i),2) == 0
      xx = [x(i):-2:2]';
      xx1 = [x(i)-1:-2:1]';
   elseif x(i) == 1
      xx = 1;
      xx1 = 1;
   else
      xx = [x(i):-2:3]';
      xx1 = [x(i)-1:-2:2]';
   end
   y(i) = prod(xx./xx1);
end
return
