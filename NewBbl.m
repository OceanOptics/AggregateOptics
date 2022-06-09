function sbb=NewBbl(nc0,FE,a,b)
% Calculate the backscattering cross section using an analytic formula
% This one is different from NewBb in that this one will use the
% precalculated coefficient by NewBbBak
% Revised on May 7, 1999. 
% Ref. Chylet, Petr
% Journal of the Optical Society of America, Vol. 63, No. 11, 1973
% Xiaodong Zhang
% May 5, 1999

N = [1:nc0]';
K = [2:2:nc0]';
L = [1:2:nc0]';
H = [1:nc0/2]';

EL1 = 2*L+1;
EL3 = L.*(L+1);
EL2 = EL1./EL3;

EK1 = 2*K+1;
EK3 = K.*(K+1);
EK2 = EK1./EK3;

aL = a(L);
bL = b(L);
FEL = FE(L);

aK = a(K);
bK = b(K);
FEK = FE(K);

ML = (-1).^H;
MK = -ML;

T1 = (2*N+1)'*(abs(a).^2+abs(b).^2);

tem = MK.*EK1./FEK;
temKa = tem.*aK;
temKb = tem.*bK;

tem = ML.*EL1.*FEL;
temLa = tem.*aL;
temLb = tem.*bL;

T2 = K;
for i=1:length(K)
   tem = (EK3(i)-EL3)';
   T2(i) = real(conj(temKa(i))./tem*temLa+conj(temKb(i))./tem*temLb);
end

clear temKa temKb temLa temLb

temK = MK.*FEL.*EL2;
temL = ML.*FEL.*EL2;
T3 = real(sum(temK.*aL)*sum(temL.*conj(bL)));

sbb = T1/2 + sum(T2) + T3;

return