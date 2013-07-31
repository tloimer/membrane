function ms = mstackstruct(theta,mem,f) %-------------------------- mstackstruct
%MSTACKSTRUCT Construct a membrane stack structure MS.
%  MSTACKSTRUCT(THETA,SUBSTANCE,MEMBRANE,FMODEL) constructs a structure which
%  describes a stack of individual membranes, each of which may consist of
%  different layers. The structure of the membrane stack is given by passing a
%  cell vector MEMBRANE of further cell vectors. The first level corresponds to
%  the separate membranes in the stack, the second level corresponds to the
%  layers in each membrane. The membrane stack structure MS is constructed in a
%  way to also accomodate the solution.
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
%    MS.colors
%    MS.printsetup
%    MS.writeflowsetups
%    MS.mfluxliquid
%    MS.freesetup
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
%  See also ASYM, FMODEL, MEMBRANE, MNUM, SUBSTANCE.

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
  'writeflowsetups',@writeflowsetups,'mfluxliquid',@mfluxliquid,...
  'freesetup',[],'substance',[],'membrane',membranes);

end %%% END MSTACKSTRUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MSTACKSTRUCT %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function [memcell, nmembranes, nlayers] = expand2cell(mem) %-------- expand2cell
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
%WRITEFLOWSETUPS Write flow setup structs to the membrane struct.
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

% With m = (kappa_i/nu) * (p_i - p_(i-1)) / L_i,   p0 |XX| p1 |XXX| p2 ... |X| pn
%						       L1,kap1  L2     .... Ln
%  m*nu* Sum_i=1^n (L_i/kappa_i) = Sum_i=1^n (p_i - p_i-1) = p_n - p_0.
sumL_kappa = 0;
for i = 1:length(ms.membrane)
  sumL_kappa = sumL_kappa + sum([ms.membrane(i).layer(:).matrix.L] ...
    ./[ms.membrane(i).layer(:).matrix.kappa]);
end
m = (p1-p2)/s.nul(T1)/sumL_kappa; % = (p1-p2)/(s.nul(T1)*sumL_kappa)
end %----------------------------------------------------------- end mfluxliquid

function printsetup(ms) %-------------------------------------------- printsetup
fprintf('  -1-\n');
nmembranes = length(ms.membrane);
for i = 1:nmembranes
  nlayers = length(ms.membrane(i).layer);
  fprintf('______ Membrane %d____________________________________________________________\n',i);
  for j = 1:nlayers
    layer = ms.membrane(i).layer(j);
    fprintf(['    Layer %d: dia = %.0f nm, thickness = %.2g mm, '...
	    '2ph-model: %s, theta = %.0f°.\n'], j, layer.matrix.dia*1e9,...
	    layer.matrix.L*1e3, layer.fmodel.name, layer.theta);
    if j < nlayers
      fprintf('    -----------------------------------------------------------------------\n');
    end
  end
  fprintf('_____________________________________________________________________________\n\n');
end
fprintf('  -2-\n');
end %------------------------------------------------------------ end printsetup
