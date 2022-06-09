

% Default line and marker size
set(0,'DefaultLineLineWidth',1.0)
set(0,'DefaultLineMarkerSize',8)
set(0,'DefaultAxesLineWidth',1.0)
set(0,'DefaultAxesFontSize',24)
set(0,'DefaultAxesFontName','Calibri')

% % Default line and marker size
% set(0,'DefaultLineLineWidth',1.25)
% set(0,'DefaultLineMarkerSize',10)
% set(0,'DefaultAxesLineWidth',1.0)
% set(0,'DefaultAxesFontSize',24)
% set(0,'DefaultAxesFontName','helvetica')


clear all
%close all

m = [1.05 1.17]
mp = [0.01 0.0001]
dens = [1.380 2.650]
ln_style = {'-','--'}
lambda = 0.660*0.75 % wavelength of incident light [micron]

D_A = logspace(log10(0.2),log10(200),400); % particle diameter [micron]
V_sphere = (pi/6)*D_A.^3;
[F, V_M] = KhelifaHill(D_A);



Cext_aggr = nan(length(m), length(D_A));
Cext_sphere = nan(length(m), length(D_A));

for km = 1:length(m)
    
    for kd = 1:length(D_A)
        kd
        [Cext_aggr(km,kd), tmp, tmp] = LatimerAggregate(D_A(kd), lambda, m(km),mp(km), F(kd));
    end

    Cext_sphere(km,:) = Cext_AD(D_A, lambda, m(km),mp(km));

end

% **Need to trace numerical reason for getting complex result from
% LatimerAggregate.m
Cext_aggr = abs(Cext_aggr);
Cext_aggr(~isfinite(Cext_aggr)) = NaN;


figure
set(gcf, 'PaperUnits','inches', 'PaperSize',[9 6], 'PaperPosition',[0.25 0.25 8.5 5.5])
                              % 'PaperSize',[W H], 'PaperPosition',[L   B   W  H]
hold on
for km = 1:length(m)
    plot(D_A, Cext_aggr(km,:)./(V_M*dens(km)), ['b' ln_style{km}])
    plot(D_A, Cext_sphere(km,:)./(V_sphere*dens(km)), ['r' ln_style{km}])
end
% set(gca, 'xscale','log','yscale','log')
% set(gca, 'xlim',[0.7 1e3], 'xtick',[1 1e1 1e2 1e3])
% set(gca, 'ylim',[1e-3 2e0], 'ytick',[1e-3 1e-2 1e-1 1])
set(gca, 'xscale','log','yscale','linear')
set(gca, 'xlim',[0.7 200], 'xtick',[1 1e1 1e2])
set(gca, 'ylim',[0 1.25], 'ytick',[0 0.25 0.5 0.75 1.0])
xlabel('D_A [\mum]')
ylabel('c_p / Mass [m^{2}/g]')

%print -depsc2 'Aggr_cpstar_vs_size.eps'


xi = 2.5:0.05:5.0

cpstar_aggr = nan(length(m), length(xi));
cpstar_sphere = nan(length(m), length(xi));



delta_DA = diff(D_A);
delta_DA = [delta_DA delta_DA(end)];
    
for km = 1:length(m)

    for kx = 1:length(xi)

        %N_D = D_A.^-xi(kx);
        % adapted from EB's Kahlifa_populations.m
%     	N_D = 4.83*1e4*D_A.^(-xi(kx));
% 		N_D = N_D*1e-12; %in order to have #particles in micron^3
%         N_D = N_D.*delta_DA;
        
        N_D = D_A.^-xi(kx) .* delta_DA;
        
        
        cpstar_aggr(km,kx) = nansum(Cext_aggr(km,:).* N_D)./nansum(V_M*dens(km).*N_D)
        cpstar_sphere(km,kx) = nansum(Cext_sphere(km,:).* N_D)./nansum(V_sphere*dens(km).* N_D)


    end

end


figure
set(gcf, 'PaperUnits','inches', 'PaperSize',[9 6], 'PaperPosition',[0.25 0.25 8.5 5.5])
                              % 'PaperSize',[W H], 'PaperPosition',[L   B   W  H]
hold on
for km = 1:length(m)
	plot(xi, cpstar_aggr(km,:), ['b' ln_style{km}])
    plot(xi, cpstar_sphere(km,:), ['r' ln_style{km}])
end
set(gca, 'xscale','linear','yscale','linear')
set(gca, 'xlim',[2.5 5.0], 'xtick',[2.5 3 3.5 4 4.5 5])
set(gca, 'ylim',[0 0.8], 'ytick',[0 0.2 0.4 0.6 0.8])
xlabel('Power law slope, J')
ylabel('c_p / Mass [m^{2}/g]')

print -depsc2 'Symp09_Aggr_cpstar_vs_PSDslope.eps'
