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
%  individual flow regimes in or in front of a layer.
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
  'printsolution',@printsolution,'singlemstofl',@singlemstofl,...
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
%  MSANY.PRINTSOLUTION(MS) Print the solution stored in MS.

% Print the substance and upstream condition
psat = ms.substance.ps(ms.T1);
fprintf(['%s, T1 = %.2f K, p1 = %.3g kPa (psat = %.3g kPa, pred = %.2f). '...
	 'm = %.3g kg/m2s.\n'],...
	ms.substance.name, ms.T1, ms.p1in/1e3, psat/1e3, ms.p1in/psat, ms.m);
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
  [isup,Tup,pup,upcolor,isfilm] = upstreamflow(ms.membrane(i));
  if isup
    printstate(Tup,pup,upcolor);
    if isfilm
      fprintf('    at liquid film, T =%+6.2f K\n',...
	      ms.membrane(i).flow(1).T(end) - ms.T2);
    end
  end
  if isup
    printstartstate(ms.membrane(i).flow(1));
  end
  fprintf('______ Membrane %d____________________________________________________________\n',i);
  nlayers = length(ms.membrane(i).layer);
  for j = 1:nlayers
    layer = ms.membrane(i).layer(j);
    printendstate(layer.flow(end));
    fprintf(['    Layer %d: dia = %.0f nm, L = %.2g mm, km = %.3g W/mK, '...
	     'theta = %.0f°.\n'], j, layer.matrix.dia*1e9,...
	    layer.matrix.L*1e3, layer.matrix.km, layer.theta);
    printstartstate(layer.flow(1));
    if j < nlayers
      fprintf('    -----------------------------------------------------------------------\n');
    end
  end
  fprintf('_____________________________________________________________________________\n\n');
end
psat = ms.substance.ps(ms.T2); % Overwrite psat
fprintf('Downstream state: T2 = %.2f K, p2 = %.3g kPa (psat = %.3g kPa, pred = %.2f)\n',...
	ms.T2, ms.p2/1e3, psat/1e3, ms.p2/psat);

%--- nested functions ------------------------------------- nested functions ---

function [isup,Tup,pup,upcolor,isfilm] = upstreamflow(amembrane) %--------------
  isfilm = false;
  if isstruct(amembrane.flow) && ~isempty(amembrane.flow(end).T)
    isup = true;
    Tup = amembrane.flow(end).T(end);
    % flow.p is silently assumed to also exist
    pup = amembrane.flow(end).p(end);
    upcolor = amembrane.flow(end).color;
    if length(amembrane.flow) == 2
      isfilm = true;
    end
  else
    isup = false;
    Tup = amembrane.layer(1).flow(end).T(end);
    pup = amembrane.layer(1).flow(end).p(end);
    upcolor = amembrane.layer(1).flow(end).color;
  end
end % --------------------------------------------------------------------------

function printendstate(flow) %--------------------------------------------------
% PRINTENDSTATE Print upstream state within a given flow region.

  % Could also do printstate(flow(end).T(end)) and call printendstate(...flow)
  % instead of printendstate(...flow(end).
  printstate(flow.T(end),flow.p(end),flow.color);
end %---------------------------------------------------------------------------

function printstartstate(flow) %------------------------------------------------
% PRINTSTARTSTATE Print downstream state within a given flow region.
  printstate(flow.T(1),flow.p(1),flow.color);
end %---------------------------------------------------------------------------

function printstate(T,p,color) %------------------------------------------------
  switch color
    case 'r'
      phasestring = 'gaseous';
    case 'b'
      phasestring = 'liquid';
    case 'g'
      phasestring = 'two-phase';
    otherwise
      phasestring = 'phase?';
  end
  fprintf('%20sT =%+6.2f K, p = %.3g kPa, pred = %.2f, %s\n', ' ',...
	  T-ms.T2, p/1e3, p/ms.substance.ps(T), phasestring);
end %---------------------------------------------------------------------------
end %--------------------------------------------------------- end printsolution
