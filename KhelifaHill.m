% ========================================================================
% Khelifa and Hill 2006 model for size-dependent fractal dimension 
%
% F = KhelifaHill(D) assumes common values for Khelifa and
% Hill model parameters. D is aggregate size in micron. F is the solid
% fraction of the aggregate.
%
% F = KhelifaHill(D, d3c,Dc,Dp) specifies the Khelifa and Hill
% model parameters. In this case, d3c is the low value of the fractal
% dimension, occuring at aggregate size Dc; Dc has units of micron; and Dp
% is the primary particle diameter in micron.
%
% Also callable as [F,Vsolid] = KhelifaHill(...) where Vsolid is the volume
% of solid mass in the aggregate in micron^3.
%
% WH Slade (wayne.slade@gmail.com) 
% ========================================================================

function [F, V_M] = KhelifaHill(D_A, varargin)

    if nargin==1
        d3c = 2     % lowest value of fractal dimension, occurs at size Dc     
        Dc = 2000	% [micron]
        Dp = 1      % primary particle diameter [micron]
    elseif nargin==4
        d3c = varargin{1};
        Dc = varargin{2};
        Dp = varargin{3};    
    else
        error('Incorrect number of function arguments.')
    end
    
    % Khelifa and Hill 2006
    beta = log10(d3c./3)./log10(Dc./Dp);
    d3 = 3.*(D_A./Dp).^beta;	% fractal dimension as function of aggregate size
    
    F = (D_A./Dp).^(d3-3);  % solid fraction as function of aggregate size

    V_M = F.*(pi/6).*D_A.^3;

    
return
    
