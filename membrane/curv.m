function curv=curv()
%CURV       Curvature [1/m].
%
%  Calls COSTHETA, EPSILON, KAPPA.

% parallel plate: kap = eps*b^2/12, b .. distance
%   curv = costheta*2/b.
% capillary bundle: kap = eps*R^2/24,
%   curv = costheta*2/R.
%
% epsilon: Lückengrad - volume fraction of pores

% capillary bundle
curv = costheta/sqrt(6*kappa/epsilon);
% parallel plates
%curv = costheta/sqrt(3*kappa/epsilon);

% HERLEITUNG R = R(kappa) für parallele Kapillaren
% siehe auch Scheidegger (1974, S. 127 ff).
% mittlere Geschwindigkeit, Hagen-Poiseuille:  u = (-dp/dx)R^2/8 mu
% superficial velocity U = u*epsilon
% Darcysches Gesetz: rho U = (-dp/dx)kappa/nu
% Identifikation => kappa = epsilon*R^2/8
% Kapillaren in alle drei Richtungen => kappa = epsilon*R^2/24
