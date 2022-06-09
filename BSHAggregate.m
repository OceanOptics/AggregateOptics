
% ========================================================================
% Boss, Slade, and Hill Aggregate Model
% Latimer aggregate model meets Khelifa and Hill
% 
% Cext = BSHAggregate(D, lambda, m,mp) assumes common values for Khelifa and
% Hill model parameters linking fractal dimension and aggregate size. D is
% aggregate size in micron; lambda is wavelength of incident light in the 
% medium in micron; and m,mp are the real and imaginary parts of the 
% relative index of refraction.
%
% Cext = BSHAggregate(D, lambda, m,mp, d3c,Dc,Dp) specifies the Khelifa and Hill
% model parameters. For this model, d3c is the low value of the fractal
% dimension, occuring at aggregate size Dc; Dc has units of micron; and Dp
% is the primary particle diameter in micron.
%
%
% This function is not vectorized. All arguments must be scalar.
%
% WH Slade (wayne.slade@gmail.com) 
% ========================================================================

function Cext_aggr = BSHAggregate(D, lambda, m,mp, varargin)

    % Callable as 4 or 7 argument function; short version assumes values 
    % for the Khelifa and Hill model parameters
    if nargin==4
        d3c = 2;     % lowest value of fractal dimension, occurs at size Dc     
        Dc = 2000;	% [micron]
        Dp = 1;      % primary particle diameter [micron]
    elseif nargin==7
        d3c = varargin{1};
        Dc = varargin{2};
        Dp = varargin{3};    
    else
        error('Incorrect number of function arguments.')
    end
    
    % Khelifa and Hill 2006
    [F, V_M] = KhelifaHill(D, d3c,Dc,Dp);

    [Cext_aggr, Cext_spheroid, Cext_shell] = LatimerAggregate(D, lambda, m,mp, F); 
    
return
    


% ========================================================================
% Derives Cext for Latimer-style aggregate model, using Shepelovich model
% for spheroid, and Zhang code for spherical shell
%   D is diameter of sphere enclosing aggregate [micron]
%   lambda is wavelength of incident light [micron]
%   m,mp are real,imaj index of refraction of constituent particles
%   F is solids fraction of aggregate
%
% Drats, this function is not vectorized. All arguments must be scalar.
%
% WH Slade (wayne.slade@gmail.com) 
% ========================================================================
function [Cext_aggr, Cext_spheroid, Cext_shell] = LatimerAggregate(D, lambda, m,mp, F)

    % (1) Spheroid model, refractive index is adjusted for solid fraction
    epsilon = 1/3;
    m_eff = F*(m-1) + 1; % Maxwell Garnet model, see Sorensen 2001 Eq 63
    mp_eff = F*mp;    
    Cext_spheroid = shepelovich_spheroid(D,lambda,m_eff,mp_eff,epsilon);

    % (2) Solid shell model
    % Outer size parameter for shell is enclosing diameter of aggregate
    y = (2*pi/lambda).*(D/2);
    % Inner size parameter, calculate based on volume of water in aggregate
    Vw = (1-F).*(4/3).*pi.*(D./2).^3;
    Dw = (6*Vw./pi).^(1/3);
    x = (2*pi/lambda).*(Dw/2);
    % Zhang's code for coated sphere
    [beta Qsca Qbb Qext] = VSFv2(1.0,x,1.15,y,[0:10:180]);
    Cext_shell = Qext.*pi.*(D./2).^2;

    % Latimer result is average of the spheroid and shell
    Cext_aggr = (Cext_spheroid+Cext_shell)/2;
    
return



% ========================================================================
% Shepelovich spheroid
% ========================================================================
function Cext = shepelovich_spheroid(D,lambda,m,mp,epsilon)

    % Size parameter of spheroid
    x = (2*pi/lambda).*(D/2);    % size parameter

    % wavelength to wavenumber
    k = 2*pi/lambda;
    
    % Integral from Shepelovich et al 2001, Eq (18)
	Cext = quadgk(@(xs)Cext_fxs(xs,x,epsilon,m,mp,k),x,x.*epsilon);

return

% For quad, implements integrated function from Shepelovich Eq (18),
% integrand xs is the size parameter of the sphere
% x is size parameter for spheroid, x = kb, k is wavenumber
% epsilon is axial ratio of spheroid, epsilon = a/b
% m,mp are real and imaginary index of refraction
% k is wavenumber [1/micron]
function g = Cext_fxs(xs, x, epsilon, m,mp, k)
    
    % Extinction of sphere acoring to anomalous diffraction model of VDH
    rho = 2*xs*(m-1);
    beta = atan(mp/(m-1));
    Qext = 2-4*exp(-rho.*tan(beta)).*(cos(beta).*sin(rho-beta)./rho+(cos(beta)./rho).^2.*cos(rho-2*beta))+4*(cos(beta)./rho).^2.*cos(2*beta);
    
    % power-law size distribution of equivalent sphere polydispersion
    f_xs = epsilon^4.*x.^5./(epsilon.^2-1)./xs.^5.*sqrt((epsilon.^2-1)./(x^2*epsilon.^2-xs.^2));
    
    g = (Qext.*pi.*x.^2./k.^2) .* f_xs;
    
return




