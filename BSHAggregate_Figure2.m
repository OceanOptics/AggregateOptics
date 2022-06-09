
clear all

% Similar to Figure 2 in "The beam attenuation of aquatic aggregates:
% implications to mass-specific attenuation" manuscript
F = [0.01 0.1 0.3 0.5 0.99]
ln_color = {'b','g','r','c','m'}
m = [1.15 1.05]
mp = [0.000 0.0001]
ln_style = {'-','--'}

lambda = 0.660*0.75 % wavelength of incident light [micron]

CEXT_VM = cell(length(m),length(F))
CEXT_sp_VM = cell(length(m),length(F));
CEXT_ho_VM = cell(length(m),length(F));

h = nan(size(F));

D_A = logspace(log10(0.7),log10(1000),50); % particle diameter[micron]

for kf = 1:length(F)
    
    V_M = F(kf).*(pi/6).*D_A.^3;
    
    for km = 1:length(mp)
        
        Cext = nan(size(D_A));
        Cext_sp = nan(size(D_A));
        Cext_ho = nan(size(D_A));
        
        
        for kd = 1:length(D_A)
            [Cext(kd), Cext_sp(kd), Cext_ho(kd)] = LatimerAggregate(D_A(kd), lambda, m(km),mp(km), F(kf));
        end
        
        CEXT_VM{km,kf} = Cext./V_M
        CEXT_sp_VM{km,kf} = Cext_sp./V_M;
        CEXT_ho_VM{km,kf} = Cext_ho./V_M;
        
    end

end



% Go back through cell arrays of results and prepare the plot
figure, hold on
for kf = 1:length(F)
    for km = length(m):(-1):1
       h(kf) = plot(D_A, CEXT_sp_VM{km,kf}, [ln_color{kf} ln_style{km}])
       plot(D_A, CEXT_ho_VM{km,kf}, [ln_color{kf} ln_style{km}])        
    end 
end
set(gca, 'xscale','log','yscale','log')
set(gca, 'xlim',[0.7 1e3])
legend(h,strcat('F=', strtrim(cellstr(num2str(F')))))
xlabel('D_A [\mum]')
ylabel('C_{ext} / V_M [\mum^{-1}]')
