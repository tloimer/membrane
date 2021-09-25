function ms = mstackstruct(theta,mem,f)
%MSTACKSTRUCT A stack of membranes, each consisting of several layers.
%  MSTACKSTRUCT(THETA,MEMBRANE,FMODEL) constructs a structure which
%  describes a stack of individual membranes, each of which may consist of
%  different layers. The structure of the membrane stack is given by
%  passing a cell vector MEMBRANE of further cell vectors. The first level
%  corresponds to the separate membranes in the stack, the second level
%  corresponds to the layers in each membrane. The membrane stack structure
%  MS is constructed in a way to also accomodate the solution.
%
%  Invoke MS = MS.WRITEFLOWSETUPS(T1,T2,SUBSTANCE,MS) or MS =
%  MS.WRITECURVSETUPS(CURV,T1,T2,SUBSTANCE,MS) after creating a MS. In the
%  latter case, THETA given above is meaningless.
%
%  The membrane stack struct MS contains the fields
%    MS.m               Mass flux [kg/m2s]
%    MS.T1              Upstream temperature [K]
%    MS.p1in            Upstream pressure, desired [Pa]
%    MS.p1sol           Upstream pressure, solution [Pa]
%    MS.a1              Upstream vapor fraction
%    MS.q1              Heat flux [W/m2]
%    MS.T2              Downstream temperature [K]
%    MS.p2              Downstream pressure [Pa]
%    MS.a2              Downstream vapor fraction
%    MS.q2              Downstream heat flux [W/m2]
%    MS.colors          Colors to plot liquid, gaseous and two-phase flow.
%    MS.printsetup      Print the membrane and layer structure
%    MS.printsolution   Print the solution, e.g., after mnumadiabat is called
%    MS.getsolution     Output the flattened solution values
%    MS.plotsolution    Plot temperature and pressure distributions
%    MS.plotT           Plot temperature distribution
%    MS.singlemstofl    Convert a MS-struct to an (obsolete) flowstruct
%    MS.writeflowsetups See MSTACKSTRUCT>WRITEFLOWSETUPS
%    MS.writecurvsetups See MSTACKSTRUCT>WRITECURVSETUPS
%    MS.mfluxliquid     Mass flux of the liquid through the membrane stack
%    MS.mfluxknudsen    Purely free molecular flux of the gas phas.
%    MS.mfluxviscous    Purely viscous mass flux of the gas phase
%    MS.freesetup       Flow setup for the free space between membranes
%    MS.substance
%    MS.membrane
%  MS.membrane is a struct of length length(MEMBRANE). It contains
%    MS.MEMBRANE.T1
%    MS.MEMBRANE.p1
%    MS.MEMBRANE.a1
%    MS.MEMBRANE.q1
%    MS.MEMBRANE.layer
%    MS.MEMBRANE.zscale
%    MS.MEMBRANE.flow
%  MS.MEMBRANE.layer is a struct of a length corresponding to the number of
%  layers in each membrane. MS.MEMBRANE(i).layer(j) contains
%    MS.MEMBRANE.LAYER.theta
%    MS.MEMBRANE.LAYER.matrix
%    MS.MEMBRANE.LAYER.fmodel
%    MS.MEMBRANE.LAYER.flsetup
%    MS.MEMBRANE.LAYER.calc
%    MS.MEMBRANE.LAYER.flow
%  MS.MEMBRANE.LAYER.flow is a struct, the length of which corresponds to
%  the individual flow regimes in a layer.
%
%  See also ASYM, FMODEL, MEMBRANE, SUBSTANCE, FLOWSETUP,
%           MSTACKSTRUCT>WRITEFLOWSETUPS, MSTACKSTRUCT>WRITEFLOWSETUPS,
%           MSTACKSTRUCT>MFLUXVISCOUS, MSTACKSTRUCT>MFLUXLIQUID,
%           MSTACKSTRUCT>MFLUXKNUDSEN, MSTACKSTRUCT>PRINTSETUP,
%           MSTACKSTRUCT>PRINTSOLUTION, MSTACKSTRUCT>PLOTSOLUTION,
%           MSTACKSTRUCT>PLOTT.

%  To plot the upstream boundary layer, for the z-coordinate use the
%  transformation (from FLOW12.m)
%    z3 = FL.flow(-FL.sol.len).z(1),  zscale = FL.sol.zscale,
%    z = z3 + zscale * log( (Fl.flow(-FL.sol.len).z-z3)/zscale + 1 ).

% expand mem to cell vector of cell vector (yes, cell vector of cell vector)
% f is expanded below
[memcell, nmembranes, nlayers] = expand2cell(mem);

% expand possible scalar respectively singleton inputs to cells of cells
if isscalar(theta)
  for i = nmembranes:-1:1
    [thetacell{i}{1:nlayers(i)}] = deal(theta);
    ntheta = nlayers;
  end
else
  [thetacell, ~, ntheta] = expand2cell(theta);
end

if length(f) == 1
  for i = nmembranes:-1:1
    [fcell{i}{1:nlayers(i)}] = deal(f);
  end
  nflay = nlayers;
else
  % expand f to cell vector of cell vector
  [fcell, ~, nflay] = expand2cell(f);
end

if ~isequal(nlayers, nflay, ntheta) % implicitly also checks nfmem, nmembranes
  error(['Membran layers %s, fmodel %s and contact angle %s arguments are '...
	 'not of the same size'], inputname(2), inputname(3), inputname(1));
end

% Allocate ms.membrane(nmembranes)
membranes(nmembranes) = struct('T1',[],'p1',[],'a1',[],'q1',[],'layer',[],...
				'zscale',[],'flow',[]);
for i = 1:nmembranes
  % Write ms.membrane(i).layer(nlayers)
  lyrs = struct('theta',thetacell{i},'matrix',memcell{i},'fmodel',fcell{i},...
		'flsetup',[],'calc',[],'flow',[]);
  membranes(i).layer = lyrs;
end

% Initialize the struct with what we know already.
ms = struct('m',[],'T1',[],'p1in',[],'p1sol',[],'a1',[],'q1',[],'T2',[],...
  'p2',[],'a2',[],'q2',[],'colors',{{'b','r','g'}},'printsetup',@printsetup,...
  'printsolution',@printsolution,'getsolution',@getsolution,...
  'plotsolution',@plotsolution,'plotT',@plotT,'plotp',@plotp,...
  'singlemstofl',@singlemstofl,'writeflowsetups',@writeflowsetups,...
  'writecurvsetups',@writecurvsetups,'mfluxliquid',@mfluxliquid,...
  'mfluxknudsen',@mfluxknudsen,'mfluxviscous',@mfluxviscous,...
  'freesetup',[],'substance',[],'membrane',membranes);

end %%% END MSTACKSTRUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MSTACKSTRUCT %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function [memcell, nmembranes, nlayers] = expand2cell(mem) %-------- expand2cell
%EXPAND2CELL Expand input to cell array of cell array.
nmembranes = length(mem);
if isstruct(mem) || isnumeric(mem)
  % nmembranes - the number of individual membranes, possibly consisting of
  %              several layers
  % nlayers -    vector of length nmembranes, number of layers, e.g. [1 3 2]
  nlayers = ones(1,nmembranes);
  for i = nmembranes:-1:1
    memcell{i} = { mem(i) };
    % now, memcell = { {mem1} {mem2} {mem3} }
  end
elseif iscell(mem)
  for i = nmembranes:-1:1
    nlayers(i) = length(mem{i});
    if isstruct(mem{i}) || isnumeric(mem{i})
      for j = nlayers(i):-1:1
	memcell{i}{j} = mem{i}(j);
      end
    elseif iscell(mem{i})
      memcell{i} = mem{i};
    end
  end
end
% Now we have a memcell cell array of structs, eg,
% { {mem11 mem12 mem13} {mem21 mem22} {mem31} }
% here, nmembranes = 3 and nlayers = [3 2 1]
end %----------------------------------------------------------- end expand2cell

function ms = writeflowsetups(T1,T2,s,ms) %--------------------- writeflowsetups
%WRITEFLOWSETUPS Add flow setup structs to a membrane stack struct.
%  MS = WRITEFLOWSETUPS(T1,T2,SUBSTANCE,MS) writes flow setup structures to
%  each layer in the membrane struct MS. The substance is written to
%  MS.substance.

% cycle through all membranes
nmembranes = length(ms.membrane);
for i = 1:nmembranes
  % and through all layers
  nlayers = length(ms.membrane(i).layer);
  for j = 1:nlayers
    % Could optimize flowsetup a bit here, by computing some stuff only once
    ms.membrane(i).layer(j).flsetup = flowsetup(T2,T1,ms.membrane(i).layer(j).theta,...
      s,ms.membrane(i).layer(j).matrix,ms.membrane(i).layer(j).fmodel);
  end
end
ms.substance = s;
ms.freesetup = flowsetup(s);
end %------------------------------------------------------- end writeflowsetups


function ms = writecurvsetups(curv,T1,T2,s,ms) %---------------- writecurvsetups
%WRITECURVSETUPS Add flow setup structs to a membrane stack struct.
%  MS = WRITECURVSETUPS(CURV,T1,T2,SUBSTANCE,MS) writes flow setup structures
%  given a radius of curvature to each layer in the membrane struct MS.
%  The substance is written to MS.substance.

% cycle through all membranes
nmembranes = length(ms.membrane);
for i = 1:nmembranes
    % and through all layers
    nlayers = length(ms.membrane(i).layer);
    for j = 1:nlayers
        % Could optimize flowsetup a bit here, by computing some stuff only once
        ms.membrane(i).layer(j).flsetup = curvsetup(T2, T1, curv, s, ...
                ms.membrane(i).layer(j).matrix, ms.membrane(i).layer(j).fmodel);
    end
end
ms.substance = s;
ms.freesetup = curvsetup(s);
end %------------------------------------------------------- end writecurvsetups

function m = mfluxliquid(T1,p1,p2,s,ms) %--------------------------- mfluxliquid
%MFLUXLIQUID Mass flux for the flow of the liquid phase.
%  MFLUXLIQUID(T1,P1,P2,SUBSTANCE,MS) returns the mass flux of the liquid
%  phase through a membrane stack MS. Does not err if the temperature T1 is
%  above the critical temperature of SUBSTANCE, but returns the mass flux
%  for the flow of the gaseous phase.

% With m = (kappa_i/nu) * (p_i - p_(i-1)) / L_i,   p0 |XX| p1 |XXX| p2 ... |X| pn
%						       L1,kap1  L2     .... Ln
%  m*nu* Sum_i=1^n (L_i/kappa_i) = Sum_i=1^n (p_i - p_i-1) = p_n - p_0.
sumL_kappa = 0;
for i = 1:length(ms.membrane)
  for j = 1:length(ms.membrane(i).layer)
    % The sum must be done manually.
    % sum([ms.membrane(i).layer(:).matrix.L] produced an error:
    %  Scalar index required for this type of multi-level indexing.
    sumL_kappa = sumL_kappa + ms.membrane(i).layer(j).matrix.L ...
		 /ms.membrane(i).layer(j).matrix.kappa;
  end
end
m = (p1-p2)/s.nul(T1)/sumL_kappa; % = (p1-p2)/(s.nul(T1)*sumL_kappa)
% No valid range checking...
% if isinf(s.ps(T1))
% a very rough estimate for the gas flux
% m = (p1-p2)/s.v(T1,(p1+p2)/2)/sumL_kappa;
end %----------------------------------------------------------- end mfluxliquid

function m = mfluxknudsen(T1,p1,p2,s,ms) %------------------------- mfluxknudsen
%MFLUXKNUDSEN Mass flux for purely molecular flow [kg/m2s].
%  MFLUXKNUDSEN(T1,P1,P2,SUBSTANCE,MS) computes the mass flux of SUBTANCE
%  for purely molecular flow through the membrane stack MS.

%   m = kappa beta Kn/nug dp/dz, Kn/nu = f(T),
%
%   m/(beta Kn/nu) * L/kappa = p1 - p2,
%
%   m = (p1 -p2) /( L_1/(beta*kappa*kn_nu)_1 + (L/beta*kappa*kn_nu)_2 + ...);
sumL_all = 0;
% kn_nu = @(T) 3*sqrt(pi/(8*s.R)) / (sqrt(T)*mem.dia);
facb_dia =  3*sqrt(pi/(8*s.R*T1));
for i = 1:length(ms.membrane)
  for j = 1:length(ms.membrane(i).layer)
    % The sum must be done manually.
    % sum([ms.membrane(i).layer(:).matrix.L] produced an error:
    %  Scalar index required for this type of multi-level indexing.
    sumL_all = sumL_all + ms.membrane(i).layer(j).matrix.L...
			  * ms.membrane(i).layer(j).matrix.dia...
			  / ms.membrane(i).layer(j).matrix.kappa / facb_dia...
			  / ms.membrane(i).layer(j).matrix.beta;
  end
end
m = (p1-p2)/sumL_all;
end %---------------------------------------------------------- end mfluxknudsen



function m = mfluxviscous(T,p1,p2,s,ms) %-------------------------- mfluxviscous
%MFLUXVISCOUS Mass flux for the viscous flow of ga.
%  MFLUXVISCOUS(T1,P1,P2,SUBSTANCE,MS) computes the mass flux for the
%  isothermal, purely viscous flow of the ideal gas, without
%  free molecular flow contribution.

% Compute a guess for the mass flux.
% With
%   m = (kappa/nu) (dp/dz),
% substituting nu = mu*v and v = RT/p, yields
%   m = (kappa*p/(mu*R*T)) dp/dz,
% integrating
%   m*L = (kappa/(mu*R*T)) (p_1^1 - p_2^2)/2.
% This is summed up over all layers.

sumL_kappa = 0;
for i = 1:length(ms.membrane)
  for j = 1:length(ms.membrane(i).layer)
    sumL_kappa = sumL_kappa + ms.membrane(i).layer(j).matrix.L ...
		 / ms.membrane(i).layer(j).matrix.kappa;
  end
end

% R*T = v*p
pm = (p1 + p2) / 2.;
m = (p1*p1 - p2*p2) / (2*s.mug(T)*s.v(T,pm)*pm*sumL_kappa);
end %---------------------------------------------------------- end mfluxviscous

function fl = singlemstofl(ms) %----------------------------------- singlemstofl
%SINGLEMSTOFL Convert a mstackstruct for a homogeneous membrane to a flowstruct.
%  FL = MSANY.SINGLEMSTOFL(MS) returns the flowstruct FL from a
%  MSTACKSTRUCT MS containing a single membrane with a single layer.
flow = [ms.membrane.layer.flow ms.membrane.flow];
fl = struct('info',struct(...
	    'T1',ms.T1,'p1',ms.p1in,'p2',ms.p2,...
	    'flsetup',ms.membrane.layer.flsetup,...
	    'membrane',ms.membrane.layer.matrix,'substance',ms.substance),...
	'sol',struct('T2',ms.T2,'len',-length(flow),'states','--'),'flow',flow);
if ~isempty(ms.membrane.zscale) && flow(end).color == 'r'
  fl.sol.zscale = ms.membrane.zscale;
  fl.sol.states = '13';
end
end %---------------------------------------------------------- end singlemstofl

function printsetup(ms) %-------------------------------------------- printsetup
%PRINTSETUP Print the membrane layer structure.
%  MSANY.PRINTSETUP(MS) prints the properties of the MSTACKSTRUCT MS.
fprintf('  -Upstream-\n\n');
nmembranes = length(ms.membrane);
for i = 1:nmembranes
  nlayers = length(ms.membrane(i).layer);
  fprintf('______ Membrane %d____________________________________________________________\n',i);
  for j = 1:nlayers
    layer = ms.membrane(i).layer(j);
    fprintf(['    Pore dia = %.0f nm, L = %.2g mm, km = %.3g W/mK, '...
	     'eps/tau = %.3g, θ = %.0f°.\n'], layer.matrix.dia*1e9,...
	    layer.matrix.L*1e3, layer.matrix.km,...
	    layer.matrix.epsilon./layer.matrix.tau, layer.theta);
    if j < nlayers
      fprintf('    -----------------------------------------------------------------------\n');
    end
  end
  fprintf('_____________________________________________________________________________\n\n');
end
fprintf('  -Downstream-\n');
end %------------------------------------------------------------ end printsetup


function printsolution(ms,pa) %----------------------------------- printsolution
%PRINTSOLUTION Print temperatures and pressures at special locations.
%  MSANY.PRINTSOLUTION(MS) Print the solution stored in the membranestruct
%  MS. The upstream and downstream states in each flow regime are printed.
%  This is equivalent to the upstream and downstream state at each front.
%
%  MSANY.PRINTSOLUTION(MS,'phase') In addition to the above, print the
%  locations of the fronts of phase change.

if nargin == 1
  pa='none';
end

% Print the substance and upstream condition
psat = ms.substance.ps(ms.T1);
fprintf(['Fluid: %s, mass flux %.4g g/m2s, p1 - p2 = %.3g kPa.\n'...
	 'Upstream state: T1 = %.2f K, p1 = %.3g kPa (psat = %.3g kPa, '...
	 'pred = %.2f).\n\n'], ms.substance.name, ms.m*1e3,...
	 (ms.p1in-ms.p2)/1e3, ms.T1, ms.p1in/1e3, psat/1e3, ms.p1in/psat);
% Print information on the iteration accuracy;
% Get the calculated solution, not the provided values.
[~,Tup,pup] = upstreamflow(ms.membrane(1));
fprintf(['Iteration accuracy: p1 - p1calc = %.2g Pa (%.2g p1),\n%20s'...
	 'T1 - T1calc = %.2g K'], ms.p1in-pup, 1-pup/ms.p1in, ' ', ms.T1-Tup);
if ms.T1 ~= ms.T2
	fprintf(' (%.2g(T1-T2)).\n\n', (ms.T1-Tup)/(ms.T1-ms.T2));
elseif ms.q2 ~= 0.
	fprintf(' (%.2g q/(m cp)).\n\n',...
		(ms.T1-Tup)*ms.m*ms.substance.cpg(ms.T2,ms.p2)/ms.q2);
else
	fprintf('.\n\n');
end

% Now loop over membranes and layers
nmembranes = length(ms.membrane);
for i = 1:nmembranes
  % Print the condition far upstream of each membrane - if there is a change with
  % respect to the values at the membrane front
  [isup,Tup,pup,aup,upcolor,isfilm] = upstreamflow(ms.membrane(i));
  if isup
    printstate(Tup,pup,aup,upcolor);
    if isfilm
      fprintf('    at liquid film, T =%+6.2f K',...
	      ms.membrane(i).flow(1).T(end) - ms.T2);
      if strcmp(pa,'phase')
        fprintf(', z = %.3g mm, z/L1 = %.2f\n',...
		ms.membrane(i).flow(1).z(end)*1e3,...
		ms.membrane(i).flow(1).z(end)/ms.membrane(i).layer(1).matrix.L);
      else
	fprintf('\n');
      end
    end
    printstartstate(ms.membrane(i).flow(1));
  end
  fprintf('______ Membrane %d_______________________________________________________\n',i);
  nlayers = length(ms.membrane(i).layer);
  for j = 1:nlayers
    layer = ms.membrane(i).layer(j);
    printendstate(layer.flow(end));
    fprintf(['    Layer %d: dia = %.0f nm, L = %.2g mm, km = %.3g W/mK, '...
	     'theta = %.0f°.\n'], j, layer.matrix.dia*1e9,...
	    layer.matrix.L*1e3, layer.matrix.km, layer.theta);
    if strcmp(pa,'phase') % a quick hack
      fprintf('      at z/L = %.2f, vapor volume fraction = %.2g\n',...
	      layer.flow(end).z(1)/layer.matrix.L, layer.flow(end).a(1));
      fprintf('      at z/L = %.2f, vapor volume fraction = %.2g\n',...
	      layer.flow(1).z(end)/layer.matrix.L, layer.flow(1).a(end));
    end
    printstartstate(layer.flow(1));
    if j < nlayers
      fprintf('    --------------------------------------------------------------------\n');
    end
  end
  fprintf('________________________________________________________________________\n\n');
end
psat = ms.substance.ps(ms.T2); % Overwrite psat
fprintf(['Downstream state: T2 = %.2f K, p2 = %.3g kPa (psat = %.3g kPa, '...
	 'pred = %.2f).\n'], ms.T2, ms.p2/1e3, psat/1e3, ms.p2/psat);

%--- nested functions ------------------------------------- nested functions ---

function printendstate(flow) %--------------------------------------------------
% PRINTENDSTATE Print upstream state within a given flow region.
  printstate(flow.T(end),flow.p(end),flow.a(end),flow.color);
end %---------------------------------------------------------------------------

function printstartstate(flow) %------------------------------------------------
% PRINTSTARTSTATE Print downstream state within a given flow region.
  printstate(flow.T(1),flow.p(1),flow.a(1),flow.color);
end %---------------------------------------------------------------------------

function printstate(T,p,a,color) %------------------------------------------------
strtwo = ' pred = %.2f, %s\n';
  switch color
    case 'r'
      phasestring = 'gaseous';
      var = p/ms.substance.ps(T);
    case 'b'
      phasestring = 'liquid';
      var = p/ms.substance.ps(T);
    case 'g'
      phasestring = 'two-phase';
      strtwo = '   a = %.3f, %s\n';
      var = a;
    otherwise
      phasestring = 'phase?';
  end
  fprintf(['%20sT =%+6.2f K, p = %.3g kPa,' strtwo], ' ',...
	  T-ms.T2, p/1e3, var, phasestring);
end %---------------------------------------------------------------------------
end %--------------------------------------------------------- end printsolution

%------------------------------------------------------------------ upstreamflow
function [isup,Tup,pup,aup,upcolor,isfilm] = upstreamflow(amembrane)
  isfilm = false;
  % isstruct(amembrane.flow) failed, because .flow may be a 0x0 struct;
  if length(amembrane.flow) > 0 && ~isempty(amembrane.flow(end).T)
    isup = true;
    Tup = amembrane.flow(end).T(end);
    % flow.p is silently assumed to also exist
    pup = amembrane.flow(end).p(end);
    aup = amembrane.flow(end).a(end);
    upcolor = amembrane.flow(end).color;
    if length(amembrane.flow) == 2
      isfilm = true;
    end
  else
    isup = false;
    Tup = amembrane.layer(1).flow(end).T(end);
    pup = amembrane.layer(1).flow(end).p(end);
    aup = amembrane.layer(1).flow(end).a(end);
    upcolor = amembrane.layer(1).flow(end).color;
  end
end % --------------------------------------------------------- end upstreamflow


function sol = getsolution(ms) %------------------------------------ getsolution
%GETSOLUTION Get solution values at fronts.
%  SOL = MSANY.GETSOLUTION(MS) gets the solution values from specific
%  locations and writes these to the struct SOL. The struct SOL contains
%  the cell arrays L, z, T, p, a, q, and Kn. Each cell in the cell arrays
%  corresponds to one membrane of an membrane stack. SOL.L refers to the
%  length of the membrane, in meters, and SOL.a is the vapor mass fraction.
%  Each cell contains a vector of property values. In essence, the FLOW
%  field in the MSTACKSTRUCT MS is flattened out and purged to only contain
%  the values at fronts of phase change and boundaries.

nmembranes = length(ms.membrane);
property = {'z', 'T', 'p', 'a', 'q', 'Kn'};
n = length(property);
% create a struct containing all empty cell arrays, equivalent to
% sol = struct('L', cell(nmembranes,1), 'z', cell(nmembranes,1),...)
sol(1).L = cell(nmembranes,1);
for l = 1:n
    sol(1).(property{l}) = cell(nmembranes,1); % initializes to empty matrices
end

% loop over all membranes
for i = 1:nmembranes
    % write the front boundary layer
    j = length(ms.membrane(i).flow);
    if j > 0 && ~isempty(ms.membrane(i).flow(end).T)
        for k = j:-1:1      % loop over the elements of the front boundary layer
            for l = 1:n     % loop over the property names
                sol.(property{l}){i} = [sol.(property{l}){i}...
                            ms.membrane(i).flow(k).(property{l})(end)...
                            ms.membrane(i).flow(k).(property{l})(1)];
            end
        end
    end

    % loop over layers
    L = 0;      % sum of layer lengths
    for j = length(ms.membrane(i).layer):-1:1
        for k = length(ms.membrane(i).layer(j).flow):-1:1
            for l = 1:n     % loop over the property names
                sol.(property{l}){i} = [sol.(property{l}){i}...
                        ms.membrane(i).layer(j).flow(k).(property{l})(end)...
                        ms.membrane(i).layer(j).flow(k).(property{l})(1)];
            end
            sol.z{i}(end-1:end) = sol.z{i}(end-1:end) + L;
        end
        l = ms.membrane(i).layer(j).matrix.L;
        sol.L{i} = [sol.L{i} l];
        L = L + l;          % sum of L - to correct the z coordinate above
    end
end

end %----------------------------------------------------------- end getsolution

function plotsolution(ms) %---------------------------------------- plotsolution
%PLOTSOLUTION Plot temperature and pressure distributions.
%  MSANY.PLOTSOLUTION(MS) plots the temperature and pressure distributions
%  for the solution stored in MS.

mark = '+';
drawingarea = get(0,'ScreenSize'); % [left bottom width height]
% left margin 57, right margin 23, top margin 23, bottom margin 27
drawingarea = drawingarea + [57 27 -80 -50];
spacing = 9;
height = floor((drawingarea(4)-spacing)/2);
bpos = drawingarea(2);
tpos = drawingarea(2) + drawingarea(4) - height;

% Count the total number of layers, including front layers.
% The number of membranes.
nmembranes = length(ms.membrane);
% This is a vector of layers in each membrane. Sum(nlayers) is the total number
% of layers.
nlayers = zeros(1,nmembranes);
isup = false(1,nmembranes);
pmin = ms.p2;
for i = 1:nmembranes
  % Here, at difference to the code in upstreamflow, the emptyness of .flow is
  % not checked; The plot-command then exits with error anyway.
  % Add the number of possible upstream layers
  if length(ms.membrane(i).flow) > 0
    nlayers(i) = nlayers(i) + 1;
    isup(i) = true;
  end
  nlayers(i) = nlayers(i) + length(ms.membrane(i).layer);
  for j = 1:length(ms.membrane(i).layer)
    pmin = min([pmin ms.membrane(i).layer(j).flow(1:end).p]);
  end
end

ntotal = sum(nlayers);
if ntotal > 2
  drawingwidth = drawingarea(3);
elseif ntotal == 2
  drawingwidth = drawingarea(3)/2;
else % ntotal == 1
  drawingwidth = drawingarea(3)/3;
end

lpos = floor(drawingarea(1) + drawingarea(3) - drawingwidth);
%drawingwidth = floor(drawingwidth);

% width of one layer = (drawingwidth - (nmembranes - 1)*spacing)/ntotal;
% width of each figure = nlayers(i)*(width of one layer);
widths = round(nlayers*(drawingwidth - (nmembranes-1)*spacing)/ntotal);

Tmax = max([ms.T1 ms.membrane(1).layer(1).flow(end).T(end)]);
Tlim = [floor(ms.T2) ceil(Tmax)];
plim = [floor(pmin/1e4) ceil(ms.p1in/1e4)]*1e4;

% Now walk through the membranes, and plot all layers
for i = 1:nmembranes
  ht = figure('Name','Temperature','OuterPosition',[lpos tpos widths(i) height]);
  hp = figure('Name','Pressure',   'OuterPosition',[lpos bpos widths(i) height]);
  lpos = lpos + widths(i) + spacing;

  % Plot the first layer, probably with the front boundary layer(s)
  jl = ms.membrane(i).layer(1);
  L = jl.matrix.L;
  if isup(i)
    first = [1 2];
    nflow = length(ms.membrane(i).flow);
    sp = subplot(1,nlayers(i),first);
    plot(sp,ms.membrane(i).flow(nflow).z/L,ms.membrane(i).flow(nflow).p/1e5,...
	 'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
    xlim(sp, [ms.membrane(i).flow(nflow).z(end)/L 1]);
    ylim(sp, plim/1e5);
    figure(ht);
    st = subplot(1,nlayers(i),first);
    plot(st, ms.membrane(i).flow(nflow).z/L,ms.membrane(i).flow(nflow).T,...
	 'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
    xlim(st, [ms.membrane(i).flow(nflow).z(end)/L 1]);
    ylim(st, Tlim);
    for k = nflow-1:-1:1
      line(st, ms.membrane(i).flow(k).z/L,ms.membrane(i).flow(k).T,...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
    end
    %figure(hp);
    %subplot(sp);
    for k = nflow-1:-1:1
      line(sp, ms.membrane(i).flow(k).z/L,ms.membrane(i).flow(k).p/1e5,...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
    end
    % The first line in the layer.
    nflow = length(jl.flow);
    line(sp, jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    %figure(ht);
    %subplot(st);
    line(st, jl.flow(nflow).z/L,jl.flow(nflow).T,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
  else
    first = 1;
    nflow = length(jl.flow);
    sp = subplot(1,nlayers(i),first);
    plot(sp, jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    ylim(sp, plim/1e5);
    xlim(sp, [0 1]);
    figure(ht);
    st = subplot(1,nlayers(i),first);
    plot(st, jl.flow(nflow).z/L,jl.flow(nflow).T,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    ylim(st, Tlim);
    xlim(st, [0 1]);
  end

  % Put the remaining flow elements in the first layer into the existing plot.
  for k = nflow-1:-1:1
    line(st, jl.flow(k).z/L,jl.flow(k).T,'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
  end
  %subplot(sp);
  for k = nflow-1:-1:1
    line(sp, jl.flow(k).z/L,jl.flow(k).p/1e5,'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
  end

  % At last, the remaining layers
  nl = length(ms.membrane(i).layer);
  if nl == 1, return; end

  for j = 2:nl
    jl = ms.membrane(i).layer(j);
    L = jl.matrix.L;
    nflow = length(jl.flow);
    figure(ht);
    st = subplot(1,nlayers(i),first(end)+j-1);
    plot(st, jl.flow(nflow).z/L,jl.flow(nflow).T,...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    ylim(st, Tlim);
    xlim(st, [0 1]);
    figure(hp);
    sp = subplot(1,nlayers(i),first(end)+j-1);
    plot(sp, jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    ylim(sp, plim/1e5);
    xlim(sp, [0 1]);
    if nflow > 1
      for k = nflow-1:-1:1
	line(sp, jl.flow(k).z/L,jl.flow(k).p/1e5,...
	     'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
      end
      %figure(ht);
      %subplot(st);
      for k = nflow-1:-1:1
	line(st, jl.flow(k).z/L,jl.flow(k).T,...
	     'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
      end
    end
  end
end

end %---------------------------------------------------------- end plotsolution

function plotT(ms,i) %---------------------------------------------------- plotT
%PLOTT      Plot temperature.
%  PLOTT(MS,I) Plot temperature distribution in the I-th membrane.

global INFO CELSIUS
if nargin == 1
  i = 1;
end

mark = 'none';

nl = length(ms.membrane(i).layer);
% the total thickness of the membrane
sumL = ms.membrane(i).layer(1).matrix.L;
if nl > 1
  offz(nl-1) = 0;
end
for j = 2:nl
  offz(j-1) = sumL;
  sumL = sumL + ms.membrane(i).layer(j).matrix.L;
end

isup = false; %false(1,nmembranes);
% Here, at difference to the code in upstreamflow, the emptyness of .flow is
% not checked; The plot-command then exits with error anyway.
% Add the number of possible upstream layers
nflow = length(ms.membrane(i).flow);
if nflow > 0
  isup = true;
  % insert 3 points between each value of z, T in the upstream boundary layer
  % the first layer might be a liquid film
  j = length(ms.membrane(i).flow(nflow).z);
  zup = interp1([0:j-1]*4, ms.membrane(i).flow(nflow).z, [0:4*j-4]);
  Tup = interp1([0:j-1]*4, ms.membrane(i).flow(nflow).T, [0:4*j-4]);
  % the zup's are negative numbers
  for k = 1:4*j-3
    if zup(1) - ms.membrane(i).zscale > zup(k)
      zup(k:4*j-3) = [];
      Tup(k:4*j-3) = [];
      break;
    end
  end
% scale back the front boundary layer
%   z3 = FL.flow(-FL.sol.len).z(1),  zscale = FL.sol.zscale,
%   z = z3 + zscale * log( (Fl.flow(-FL.sol.len).z-z3)/zscale + 1 ).
  zup = zup(1) + ms.membrane(i).zscale * log( (zup-zup(1))/ms.membrane(i).zscale + 1);
end

if ms.T1 == ms.T2 || (isscalar(CELSIUS) && CELSIUS == 1)
	mkTdim = @(T) T - 273.15;
	dimensional = 1;
else
	mkTdim = @(T) (T - ms.T2)/(ms.T1 - ms.T2);
	dimensional = 0;
end

if INFO == 1
	fprintf('T1 = %.2f K, T2 = %.2f K.\n', ms.T1, ms.T2);
	fmt = '%5.3f %g\n';
end

  % Plot the first layer, probably with the front boundary layer(s)
  ht = figure('Name','Global Temperature');
  jl = ms.membrane(i).layer(1);
  if isup
    plot(zup/sumL, mkTdim(Tup),...
	'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
    if INFO == 1
      fprintf('z/L    T/K\n');
      fprintf(fmt, [zup/sumL; Tup]);
      fprintf('\n');
    end

    xlim([zup(end)/sumL 1]);
    for k = nflow-1:-1:1
      line(ms.membrane(i).flow(k).z/sumL, mkTdim(ms.membrane(i).flow(k).T),...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
      if INFO == 1
        fprintf(fmt, ...
	    [ms.membrane(i).flow(k).z/sumL; ms.membrane(i).flow(k).T]);
	fprintf('\n');
      end
    end
    % The first line in the layer.
    nflow = length(jl.flow);
    line(jl.flow(nflow).z/sumL, mkTdim(jl.flow(nflow).T),...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    if INFO == 1
      fprintf(fmt, [jl.flow(nflow).z/sumL; jl.flow(nflow).T]);
      fprintf('\n');
    end
  else
    nflow = length(jl.flow);
    plot(jl.flow(nflow).z/sumL, mkTdim(jl.flow(nflow).T),...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    xlim([0 1]);
    if INFO == 1
      fprintf(fmt, [jl.flow(nflow).z/sumL; jl.flow(nflow).T]);
      fprintf('\n');
    end
  end

  xlabel('z/(sum L)');
  if dimensional == 0
    ylim([-0.1 1.1]);
    ylabel('(T-T_2)/(T_1-T_2)');
  else
    ylabel('T [°C]');
  end

  % Put the remaining flow elements in the first layer into the existing plot.
  for k = nflow-1:-1:1
    line(jl.flow(k).z/sumL, mkTdim(jl.flow(k).T),'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
    if INFO == 1
      fprintf(fmt, [jl.flow(k).z/sumL; jl.flow(k).T]);
      fprintf('\n');
    end
  end

  % At last, the remaining layers
  if nl == 1, return; end

  for j = 2:nl
    jl = ms.membrane(i).layer(j);
    nflow = length(jl.flow);
    line([1 1]*offz(j-1)/sumL, [0.2 0.8]);
    for k = nflow:-1:1
      line((offz(j-1) + jl.flow(k).z)/sumL,mkTdim(jl.flow(k).T),...
	   'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
    end
    if INFO == 1
      fprintf(fmt, [(offz(j-1) + jl.flow(k).z)/sumL; jl.flow(k).T]);
      fprintf('\n');
    end
  end

  % Boundaries between layers
  k = ylim;
  for j = 2:nl
    line([1 1]*offz(j-1)/sumL, k(1) + [0.2 0.8]*(k(2) - k(1)));
  end

end %----------------------------------------------------------------- end plotT

function plotp(ms,i) %---------------------------------------------------- plotp
%PLOTP      Plot pressure distribution.
%  PLOTP(MS,I) Plot pressure distribution in the I-th membrane.

global INFO
if nargin == 1
  i = 1;
end

mark = 'none';
scalep = 1e-3;		% Plot kPa, not Pa.

nl = length(ms.membrane(i).layer);
% the total thickness of the membrane
sumL = ms.membrane(i).layer(1).matrix.L;
if nl > 1
  offz(nl-1) = 0;
end
for j = 2:nl
  offz(j-1) = sumL;
  sumL = sumL + ms.membrane(i).layer(j).matrix.L;
end

isup = false; %false(1,nmembranes);
% Here, at difference to the code in upstreamflow, the emptyness of .flow is
% not checked; The plot-command then exits with error anyway.
% Add the number of possible upstream layers
nflow = length(ms.membrane(i).flow);
% %{
if nflow > 0
  isup = true;
  % add one point far upstream
  j = length(ms.membrane(i).flow(nflow).z);
  zup = ms.membrane(i).flow(nflow).z;
  pup = ms.membrane(i).flow(nflow).p;
  % remove values that would yield the log of a negative number below
  % note: the zup's are negative numbers
  for k = 1:j
    if zup(1) - ms.membrane(i).zscale > zup(k)
      zup(k:j) = [];
      break;
    end
  end
% %}
  % scale back the front boundary layer
  %   z3 = FL.flow(-FL.sol.len).z(1),  zscale = FL.sol.zscale,
  %   z = z3 + zscale * log( (Fl.flow(-FL.sol.len).z-z3)/zscale + 1 ).
  zup = zup(1) + ms.membrane(i).zscale * ...
	log( (zup-zup(1))/ms.membrane(i).zscale + 1);
  % add one point far upstream
  zup = [zup -sumL];
  pup = [pup ms.membrane(i).p1];
end

ht = figure('Name','Pressure distribution');
axes(ht, 'Box', 'on');
jl = ms.membrane(i).layer(1);
if isup
  line(zup/sumL, pup*scalep,...
	'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
  if INFO == 1
    fmt = '%5.3f %g\n';
    fprintf('z/L    p/Pa\n');
    fprintf(fmt, [zup/sumL; pup]);
    fprintf('\n');
  end

  xlim([zup(end-1)/sumL 1]);
  for k = nflow-1:-1:1
    line(ms.membrane(i).flow(k).z/sumL, ms.membrane(i).flow(k).p*scalep,...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
    if INFO == 1
      fprintf(fmt, ...
	    [ms.membrane(i).flow(k).z/sumL; ms.membrane(i).flow(k).p]);
	fprintf('\n');
    end
  end
else
  xlim([0 1]);
end
nflow = length(jl.flow);
%ylim([-0.1 1.1]);
xlabel('z/L');
ylabel('p [kPa]');

% The distribution in  the first layer.
for k = nflow:-1:1
  line(jl.flow(k).z/sumL, jl.flow(k).p*scalep, 'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
  if INFO == 1
    fprintf(fmt, [jl.flow(k).z/sumL; jl.flow(k).p]);
    fprintf('\n');
  end
end

% At last, the remaining layers
if nl == 1, return; end

for j = 2:nl
  jl = ms.membrane(i).layer(j);
  nflow = length(jl.flow);
  for k = nflow:-1:1
    line((offz(j-1) + jl.flow(k).z)/sumL, jl.flow(k).p*scalep,...
	   'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
  end
  if INFO == 1
    fprintf(fmt, [(offz(j-1) + jl.flow(k).z)/sumL; jl.flow(k).p]);
    fprintf('\n');
  end
end

end %----------------------------------------------------------------- end plotp
