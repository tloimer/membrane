function pTplot(fs)
%PTPLOT     Plot a p-T Diagram.
%  PTPLOT(FLOWSTRUCT) plots the path of the process, the path of the
%  linear theory and the vapor pressure over the temperature.
%
%  See also FLOWSTRUCT.

T0 = fs.info.T0;
p0 = fs.info.p0;

% the region of interest
dT = fs.info.T0 - fs.sol.Te;
Tmin = fs.sol.Te - dT/2;
Tmax = fs.info.T0 + dT/2;
Trange = Tmin : dT/15 : Tmax;

% the complete path; add the initial point of saturated vapor
T = [fs.flow(:).T fs.info.T0];
p = [fs.flow(:).p fs.info.p0];

% the first point at the front (after condensation)
Tf = T(end-1);
pf = p(end-1);

% linear path
if fs.lin.x3==0
  Tlin = [T0 fs.lin.T4 fs.lin.Te fs.lin.Te];
  plin = [p0 p0 fs.lin.p6 p0-fs.info.dp];
else
  Tlin = [T0 fs.lin.Te fs.lin.Te];
  plin = [p0 fs.lin.p6 p0-fs.info.dp];
end

% the vapor pressure line
prange = ps(Trange);
plot(Trange,prange,'LineWidth',2,'Color',[0.8 0.8 0.8]);

% nonlinear path
line(T,p,'Color','k');

% linear path
line(Tlin,plin,'Color','k','LineStyle',':','Marker','x');
