function m = membrane(dia,epsilon,km,tname,tau,beta,L)
%MEMBRANE   Membrane properties.
%  M = MEMBRANE(DIA,EPSILON,KM,TNAME,TAU,BETA,L) returns a struct that contains
%  the membrane properties, depending on the pore diameter DIA, the void
%  fraction EPSILON, the thermal conductivity of the membrane material KM, the
%  topology of the pore space TNAME, tortuosity TAU, molecular flow correction
%  factor BETA and membrane thickness L. The topology TNAME can be one of
%  'porousround', 'porousslit', 'tube' or 'channel'.
%
%  M contains the fields
%    M.tname        Topology.
%    M.dia          Pore diameter [m].
%    M.epsilon      Void fraction, porosity.
%    M.km           Thermal conductivity of the membrane [W/mK].
%    M.tau          Tortuosity.
%    M.beta         Molecular flow factor.
%    M.L            Thickness of the membrane [m].
%    M.kappa        Permeability [m2].
%
%  Functions
%    M.fcurv        Curvature(costheta) [1/m].
%    M.fkappa       Permeability(dia,epsilon) [m2].
%    M.fdia         Diameter(kappa,epsilon) [m].
%
%  The functions FKAPPA and FDIA depend on the topology TNAME and on TAU while
%  FCURV additionaly depends on DIA.

% Struct constructor.
m = struct('tname',tname,'dia',dia','epsilon',epsilon,'km',km,'tau',tau, ...
  'beta',beta,'L',L,'kappa',[],'fcurv',[],'fkappa',[],'fdia',[]);

% thermal conductivity of the membrane material [W/mK];
% void fraction, tortuosity, beta and topology of the membrane
% tortuosity, beta: somehow a response to viscous and molecular flow,
% respectively.
% km = 1.38; % Vycor glass (Elmer, 1992).
% km = 23; % ceramic
% km = 0.22; beta = 4.67; % Celgard
% beta: molecular flow correction factor, eq. (3);
% beta = 4.67; % Celgard; beta = 8.1; % tube
%km = 0.22; epsilon = 0.39; tau = 1.3; beta = 8.1;

% TOPOLOGY
% HERLEITUNG R = R(kappa) für parallele Kapillaren
% siehe auch Scheidegger (1974, S. 127 ff).
% mittlere Geschwindigkeit, Hagen-Poiseuille:  u = (-dp/dx)R^2/8 mu
% superficial velocity U = u*epsilon
% Darcysches Gesetz: rho U = (-dp/dx)kappa/nu
% Identifikation => kappa = epsilon*R^2/8 (für parallele Kapillaren)
% Kapillaren in alle drei Richtungen => kappa = epsilon*R^2/24

% Kanalströmung (entspricht parallelen Schlitzen): u = (-dp/dx) D^2/12 mu
% Identifikation => kappa = epsilon*D^2/12

% Zusammenfassung, berücksichtige tau:
% kappa = (epsilon/tau) D^2/12 für Schlitze, Breite D, pcap = sig/(D/2)
% kappa = (epsilon/tau) R^2/24 für Poren, Radius R, pcap = 2 sig/R
% kappa = (epsilon/tau) R^2/8, tau = 1 für Röhren, pcap = 2 sig/R

% kappa = fackappa*epsilon*D^2;
% curv = faccurv/D;
% diameter =  sqrt(kappa/(fackappa*epsilon)

switch(tname)
case 'porousround'
fackappa = 1/(tau*96);
faccurv = 4;
case {'porousslit','channel'}
if strcmp(name,'channel'); tau =1;
end
fackappa = 1/(tau*12);
faccurv = 2;
case 'tube'
fackappa = 1/(tau*32);
faccurv = 4;
otherwise
error('Topologies must be one of porousround, porousslit, channel, tube.')
end
%end case name

m.fkappa = @(diameter,eps) fackappa*eps.*diameter.^2;
m.fcurv = @(costheta) faccurv*costheta./dia;
% I could also export fcurvdia(dia,costheta) ...
%m.fcurv = @(diameter,costheta) faccurv*costheta./diameter;
m.fdia = @(kappa,eps) sqrt(kappa./(fackappa*eps));

m.kappa = m.fkappa(dia,epsilon);