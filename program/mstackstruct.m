function ms = mstackstruct(theta,mem,f) %-------------------------- mstackstruct
%MSTACKSTRUCT Construct a membrane stack structure MS.
%  MSTACKSTRUCT(THETA,MEMBRANE,FMODEL) constructs a structure which describes a
%  stack of individual membranes, each of which may consist of different layers.
%  The structure of the membrane stack is given by passing a cell vector MEMBRANE
%  of further cell vectors. The first level corresponds to the separate membranes
%  in the stack, the second level corresponds to the layers in each membrane. The
%  membrane stack structure MS is constructed in a way to also accomodate the
%  solution.
%
%  Invoke MS = MS.WRITEFLOWSETUPS(T1,T2,SUBSTANCE,MS) after creating a MS.
%
%  The membrane stack structure MS contains the fields
%    MS.m
%    MS.T1
%    MS.p1in
%    MS.p1sol
%    MS.a1
%    MS.q1
%    MS.T2
%    MS.p2
%    MS.a2
%    MS.q2
%    MS.colors          Colors to plot liquid, gaseous and two-phase flow.
%    MS.printsetup      Print the membrane and layer structure.
%    MS.printsolution   Print the solution, e.g., after mnumadiabat is called.
%    MS.plotsolution    Plot temperature and pressure distributions
%    MS.singlemstofl    Convert a MS-struct to an (obsolete) flowstruct.
%    MS.writeflowsetups See MSTACKSTRUCT>WRITEFLOWSETUPS.
%    MS.mfluxliquid     Mass flux of the liquid through the membrane stack.
%    MS.freesetup       Flow setup for the free space between membranes.
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
%  MS.MEMBRANE.LAYER.flow is a struct, the length of which corresponds to the
%  individual flow regimes in a layer.
%
%  See also ASYM, FMODEL, MEMBRANE, MSTACKSTRUCT>WRITEFLOWSETUPS, SUBSTANCE.

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
  'printsolution',@printsolution,'plotsolution',@plotsolution,'singlemstofl',@singlemstofl,...
  'writeflowsetups',@writeflowsetups,'mfluxliquid',@mfluxliquid,...
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
%WRITEFLOWSETUPS Compute and write flow setup structs to the membrane struct.
%  MS = WRITEFLOWSETUPS(T1,T2,SUBSTANCE,MS) writes flow setup structures to each
%  layer in the membrane struct MS. The substance is written to ms.substance.

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

function m = mfluxliquid(T1,p1,p2,s,ms) %--------------------------- mfluxliquid
%MFLUXLIQUID Mass flux for the flow of liquid through the membrane stack.
%  M = MFLUXLIQUID(T1,P1,P2,SUBSTANCE,MS)

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
end %----------------------------------------------------------- end mfluxliquid

function fl = singlemstofl(ms) %----------------------------------- singlemstofl
%SINGLEMSTOFL Convert a mstackstruct for a homogeneous membrane to a flowstruct.
%  FL = MSANY.SINGLEMSTOFL(MS) returns the flowstruct FL from a MSTACKSTRUCT MS
%  containing a single membrane with a single layer.
flow = [ms.membrane.layer.flow ms.membrane.flow];
fl = struct('info',struct('membrane',ms.membrane.layer.matrix),...
  'sol',struct('len',-length(flow),'states','--'),'flow',flow);
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
    fprintf(['    Layer %d: dia = %.0f nm, L = %.2g mm, km = %.3g W/mK, '...
	     'theta = %.0f°.\n'], j, layer.matrix.dia*1e9,...
	    layer.matrix.L*1e3, layer.matrix.km, layer.theta);
    if j < nlayers
      fprintf('    -----------------------------------------------------------------------\n');
    end
  end
  fprintf('_____________________________________________________________________________\n\n');
end
fprintf('  -Downstream-\n');
end %------------------------------------------------------------ end printsetup


function printsolution(ms) %-------------------------------------- printsolution
%PRINTSOLUTION Print temperatures and pressures at special locations.
%  MSANY.PRINTSOLUTION(MS) Print the solution stored in the membranestruct MS.
%  The upstream and downstream states in each flow regime are printed. This is
%  equivalent to the upstream and downstream state at each front.

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
	 'T1 - T1calc = %.2g K (%.2g(T1-T2)).\n\n'],...
	ms.p1in-pup, 1-pup/ms.p1in, ' ', ms.T1-Tup, (ms.T1-Tup)/(ms.T1-ms.T2));

% Now loop over membranes and layers
nmembranes = length(ms.membrane);
for i = 1:nmembranes
  % Print the condition far upstream of each membrane - if there is a change with
  % respect to the values at the membrane front
  [isup,Tup,pup,aup,upcolor,isfilm] = upstreamflow(ms.membrane(i));
  if isup
    printstate(Tup,pup,aup,upcolor);
    if isfilm
      fprintf('    at liquid film, T =%+6.2f K\n',...
	      ms.membrane(i).flow(1).T(end) - ms.T2);
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

function plotsolution(ms) %---------------------------------------- plotsolution
%PLOTSOLUTION Plot temperature and pressure distribution.

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
%  pmin = min([pmin ms.membrane(i).layer(j).flow(1:end).p(1)]);
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
    sp = subplot(1,nlayers(i),first);
    nflow = length(ms.membrane(i).flow);
    plot(ms.membrane(i).flow(nflow).z/L,ms.membrane(i).flow(nflow).p/1e5,...
	 'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
    xlim([ms.membrane(i).flow(nflow).z(end)/L 1]);
    ylim(plim/1e5);
    figure(ht);
    st = subplot(1,nlayers(i),first);
    plot(ms.membrane(i).flow(nflow).z/L,ms.membrane(i).flow(nflow).T,...
	 'Color',ms.membrane(i).flow(nflow).color,'LineStyle','-','Marker',mark);
    xlim([ms.membrane(i).flow(nflow).z(end)/L 1]);
    ylim(Tlim);
    for k = nflow-1:-1:1
      line(ms.membrane(i).flow(k).z/L,ms.membrane(i).flow(k).T,...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
    end
    figure(hp);
    subplot(sp);
    for k = nflow-1:-1:1
      line(ms.membrane(i).flow(k).z/L,ms.membrane(i).flow(k).p/1e5,...
	   'Color',ms.membrane(i).flow(k).color,'LineStyle','-','Marker',mark);
    end
    % The first line in the layer.
    nflow = length(jl.flow);
    line(jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    figure(ht);
    subplot(st);
    line(jl.flow(nflow).z/L,jl.flow(nflow).T,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
  else
    first = 1;
    sp = subplot(1,nlayers(i),first);
    nflow = length(jl.flow);
    plot(jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    ylim(plim/1e5);
    xlim([0 1]);
    figure(ht);
    st = subplot(1,nlayers(i),first);
    plot(jl.flow(nflow).z/L,jl.flow(nflow).T,'Color',jl.flow(nflow).color,...
	 'LineStyle','-','Marker',mark);
    ylim(Tlim);
    xlim([0 1]);
  end

  % Put the remaining flow elements in the first layer
  for k = nflow-1:-1:1
    line(jl.flow(k).z/L,jl.flow(k).T,'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
  end
  subplot(sp);
  for k = nflow-1:-1:1
    line(jl.flow(k).z/L,jl.flow(k).p/1e5,'Color',jl.flow(k).color,...
	 'LineStyle','-','Marker',mark);
  end

  % At last, the remaining layers
  nl = length(ms.membrane(i).layer);
  if nl == 1, return; end

  for j = 2:nl
    figure(ht);
    st = subplot(1,nlayers(i),first(end)+j-1);

    jl = ms.membrane(i).layer(j);
    L = jl.matrix.L;
    nflow = length(jl.flow);
    plot(jl.flow(nflow).z/L,jl.flow(nflow).T,...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    ylim(Tlim);
    xlim([0 1]);
    figure(hp);
    sp = subplot(1,nlayers(i),first(end)+j-1);
    plot(jl.flow(nflow).z/L,jl.flow(nflow).p/1e5,...
	 'Color',jl.flow(nflow).color,'LineStyle','-','Marker',mark);
    ylim(plim/1e5);
    xlim([0 1]);
    if nflow > 1
      for k = nflow-1:-1:1
	line(jl.flow(k).z/L,jl.flow(k).p/1e5,...
	     'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
      end
      figure(ht);
      subplot(st);
      for k = nflow-1:-1:1
	line(jl.flow(k).z/L,jl.flow(k).T,...
	     'Color',jl.flow(k).color,'LineStyle','-','Marker',mark);
      end
    end
  end
end

end %---------------------------------------------------------- end plotsolution
