function state = downstreamstate(T2,p2,a2,q2,s,m) %------------- downstreamstate
%DOWNSTREAMSTATE Returns a struct STATE, probably calculating missing values.
%
%  STATE = DOWNSTREAMSTATE(T2,P2,A2,Q2,SUBSTANCE,MASSFLUX) constructs and returns
%  a struct STATE that describes the downstream state. If A2 is empty, the state
%  of the fluid is calculated by comparing P2 with the saturation pressure,
%  i.e., assuming free space. A saturated liquid or saturated vapor can be given
%  by setting A2 to 0 or 1, respectively, ommitting P2. For two-phase flow, the
%  homogeneous flow model FMODEL('plug') is used. Does not need MEMBRANE('free')
%  or FLOWSETUP(S) for the free space. (Why did i write those?)
%
%  If A2 is empty and P2 = PSAT, a state of a saturated vapor is returned.
%
%  The struct STATE also contains functions to set a state, see
%  DOWNSTREAMSTATE>STATESTRUCT.
%
%  DOWNSTREAMSTATE does not check the input.
%
%  See also SUBSTANCE, DOWNSTREAMSTATE>STATESTRUCT.

alphaliquid = 0;
alphavapor = 1;

% Check if the the fluid is a liquid or a vapor, if a2 is not given.
if isempty(a2)
  if p2 > s.ps(T2)
    % liquid phase
    a2 = alphaliquid;
  else
    % gaseous phase, vapor
    a2 = alphavapor;
  end
elseif isempty(p2)
  % for a two-phase mixture, assume free space
  p2 = s.ps(T2);
  if ~isfinite(p2), error('Downstream pressure out of range'); end
end

% Phase letters
gaseous = 'g';
liquid = 'l';
twophase = '2';

switch a2
  case alphavapor
    state = avapor(T2,p2,q2);
  case alphaliquid
    state = aliquid(T2,p2,q2);
  otherwise
    % p2 == s.ps(T2);
    homogeneous = fmodel('plug');
    % doth is relative to the enthalpy of the liquid at T2!
    hvap = s.hvap(T2);
    doth = q2 + m*homogeneous.xdot(a2,s.v(T2,p2),1/s.rho(T2))*hvap;
    state = atwophase(T2,doth,hvap,p2,0,0);
    % state.pliq is needed for the front, but not at the beginning (downstream
    % end) of the two-phase integrator. Therefore, atwophase does not set pliq.
    state.pliq = p2;
    state.q = q2; % must be zero here, though. Or, to simulate a two-dimensional
    % heat loss?
end

end %%% END DOWNSTREAMSTATE %%%%%%%%%%%%%%%%%%%%%%%%%%%% END DOWNSTREAMSTATE %%%

%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function state = avapor(T,p,q) %----------------------------------------- avapor
% AVAPOR    Construct and return the struct STATE of a vapor.
%
%  STATE = AVAPOR(T,P,Q)
  state = statestruct;
  state.T = T;
  state.p = p;
  state.a = 1;
  state.q = q;
  state.phase = 'g';
end %---------------------------------------------------------------- end avapor

function state = aliquid(T,p,q) %--------------------------------------- aliquid
% ALIQUID   Construct and return the struct STATE of a liquid.
%
%  STATE = ALIQUID(T,P,Q)
  state = statestruct;
  state.T = T;
  state.p = p;
  state.a = 0;
  state.q = q;
  state.phase = 'l';
end %--------------------------------------------------------------- end aliquid

function state = atwophase(T,doth,hvapK,pk,dpk,dpcap) %--------------- atwophase
% ATWOPHASE Construct and return the struct STATE of a two-phase mixture.
%
%  STATE = ATWOPHASE(T,DOTH,HVAPK,PK,DPK,DPCAP)
  state = statestruct;
  state.T = T;
  state.doth = doth;
  state.phase = '2';
  state.hvapK = hvapK;
  state.pk = pk;
  state.dpk = dpk;
  state.dpcap = dpcap;
end %------------------------------------------------------------- end atwophase

function state = statestruct %-------------------------------------- statestruct
%STATESTRUCT Construct a STATESTRUCT.
%
%  STATE = STATESTRUCT constructs the struct STATE with the fields
%    STATE.T        Temperature [K]
%    STATE.p        Pressure [Pa]
%    STATE.q        Heat flux [W/m2]
%    STATE.doth     Flux of enthalpy [W/m2]
%    STATE.phase    Phase letter: 'l', 'g' or '2'.
%    STATE.hvapK    Enthaply of vaporization at curved interface
%    STATE.pk       Pressure of the vapor, in a two-phase mixture [Pa]
%    STATE.pliq     Pressure of the liquid, in a two-phase mixture [Pa]
%    STATE.dpk      d pk/dT
%    STATE.dpcap    d pcap/ dT, where pcap = pvap - pliq
%    STATE.avapor   Function ro return a vapor state, AVAPOR(T,P,Q)
%    STATE.aliquid  Function ro return a liquid state, ALIQUID(T,P,Q)
%    STATE.atwophase Function ro return a two-phase state, see
%                    DOWNSTREAMSTATE>ATWOPHASE

state = ...
  struct('T',[],'p',[],'a',[],'q',[],'doth',[],'phase',[],'hvapK',[],'pk',[],...
    'pliq',[],'dpk',[],'dpcap',[],'avapor',@avapor,'aliquid',@aliquid,...
    'atwophase',@atwophase);
end %----------------------------------------------------------- end statestruct
