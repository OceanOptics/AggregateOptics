function Cext = Cext_AD(D, lambda, m,mp)

    k = 2*pi/lambda;
    x = (pi/lambda).*D;

    % Extinction of sphere according to anomalous diffraction model of VDH
    rho = 2*x*(m-1);
    beta = atan(mp/(m-1));
    Qext = 2-4*exp(-rho*tan(beta)) .* (cos(beta)*sin(rho-beta)./rho+(cos(beta)./rho).^2 .* cos(rho-2*beta)) ...
        + 4*(cos(beta)./rho).^2 * cos(2*beta);
    
    Cext = Qext.*(pi/4).*D.^2;
