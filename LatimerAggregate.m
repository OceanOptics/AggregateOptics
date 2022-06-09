
% ========================================================================
% Derives Cext for Latimer-style aggregate model, using Shepelovich model
% for spheroid, and Zhang code for spherical shell
%   D_A is diameter of sphere enclosing aggregate [micron]
%   lambda is wavelength of incident light [micron]
%   m,mp are real,imaj index of refraction of constituent particles
%   F is solids fraction of aggregate
%
% This function is not vectorized, all arguments must be scalar.
%
% WH Slade (wayne.slade@gmail.com) 
% ========================================================================
function [Cext, Cext_spheroid, Cext_hollow] = LatimerAggregate(D_A, lambda, m,mp, F)

    % HOLLOW SPHERE MODEL
    % Outer size parameter for shell is enclosing diameter of aggregate
    y = (pi/lambda)*(D_A);
    % Inner size parameter, calculate based on volume of water in aggregate
    % Vcore = (1-F)*(pi/6)*D_A^3
    x = (pi/lambda)*(1-F)^(1/3)*D_A;
    % Use Zhang's code for coated sphere
    [beta Qsca Qbb Qext] = VSFv2(1,x,m+sqrt(-1)*mp,y);
    Cext_hollow = Qext * (pi/4)*D_A^2;

    
    % SPHEROID MODEL
    % Latimer (1985): The gross volume of the homogeneous spheroid
    % is also set equal to the gross volume of the aggregate. Its
    % refractive index is the average value weighted in proportion
    % to the volumes of the component particles and
    % the spaces. For aggregates of three or more particles
    % modeled as prolate spheroids, the axial ratio was taken
    % as 3:1.    
    epsilon = 1/3;
    m_eff = F*(m-1) + 1; % Simplification of Gladstone-Dale rule, see see Jonasz and Fournier 2007 p466, see Sorensen 2001 Eq 63
    mp_eff = F*mp;    
    
    % "Gross volume of spheroid" should equal "gross volume" of aggregate.
    % I interpret this as the spheroid having volume equal to that of the
    % sphere enclosing aggregate. Determine the size parameter alpha=kb
    % which gives spheroid with volume equal to volume of
    % aggregate-enclosing sphere.  
    k = 2*pi/lambda;
    alpha = (k/2)*D_A*epsilon^(-2/3); % Size parameter of equivalent-volume spheroid

    % Integral from Shepelovich et al 2001, Eq (18)
	Cext_spheroid = quad(@(x)Cext_fx(x,alpha,epsilon,m_eff,mp_eff,k),alpha,alpha*epsilon);
	%Cext_spheroid = quadgk(@(x)Cext_fx(x,alpha,epsilon,m_eff,mp_eff,k),alpha,alpha*epsilon);

 
    % Latimer result is average of the spheroid and hollow sphere
    Cext = (Cext_spheroid+Cext_hollow)/2;
    
return


% For quad, implements integrated function from Shepelovich Eq (18),
% integrand x is the size parameter of the sphere
% alpha is size parameter for spheroid, alpha = kb, k is wavenumber
% epsilon is axial ratio of spheroid, epsilon = a/b
% m,mp are real and imaginary index of refraction
% k is wavenumber [1/micron]
function g = Cext_fx(x, alpha, epsilon, m,mp, k)
    
    % Extinction of sphere according to anomalous diffraction model of VDH
    rho = 2*x*(m-1);
    beta = atan(mp/(m-1));
    Qext = 2-4*exp(-rho*tan(beta)) .* (cos(beta)*sin(rho-beta)./rho+(cos(beta)./rho).^2 .* cos(rho-2*beta)) ...
        + 4*(cos(beta)./rho).^2 * cos(2*beta);
    
    % power-law size distribution of equivalent sphere polydispersion
    f_x = (epsilon^4*alpha^5)./((epsilon^2-1)*x.^5) .* sqrt((epsilon^2-1)./(alpha^2*epsilon^2-x.^2));
    
    g = (Qext.*pi.*x.^2./k^2) .* f_x;
    
return




