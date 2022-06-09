function nc=Nstop(x)
% Calculate the total steps needed for Bessel function to converge
x=abs(x);
nc = ceil(x+4.05*(x^(1/3))+2);
if rem(nc,2) ~= 0
   nc = nc+1;
end
