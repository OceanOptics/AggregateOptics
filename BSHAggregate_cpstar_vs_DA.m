

% Default line and marker size
set(0,'DefaultLineLineWidth',1.0)
set(0,'DefaultLineMarkerSize',8)
set(0,'DefaultAxesLineWidth',1.0)
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultAxesFontName','Times New Roman')


clear all
%close all

m = [1.05 1.17]
mp = [0.01 0.0001]
ln_style = {'-','--'}
lambda = 0.660*0.75 % wavelength of incident light [micron]

D_A = logspace(log10(0.7),log10(1000),80); % particle diameter [micron]

[F, V_M] = KhelifaHill(D_A);



figure
set(gcf, 'PaperUnits','inches', 'PaperSize',[5 4.0], 'PaperPosition',[0.25 0.25 4.5 3.5])
                              % 'PaperSize',[W H], 'PaperPosition',[L   B     W  H]
hold on

for km = 1:length(m)
    
    Cext_aggr = nan(size(D_A));
    for kd = 1:length(D_A)
        kd
        [Cext_aggr(kd), tmp, tmp] = LatimerAggregate(D_A(kd), lambda, m(km),mp(km), F(kd));
    end

    Cext_sphere = Cext_AD(D_A, lambda, m(km),mp(km));
    V_sphere = (pi/6)*D_A.^3;

    plot(D_A, Cext_aggr./V_M, ['b' ln_style{km}])
    plot(D_A, Cext_sphere./V_sphere, ['r' ln_style{km}])

end

set(gca, 'xscale','log','yscale','log')
set(gca, 'xlim',[0.7 1e3])
xlabel('D_A [\mum]')
ylabel('C_{ext} / V_M [\mum^{-1}]')

