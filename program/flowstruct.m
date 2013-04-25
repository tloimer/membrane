function fl = flowstruct(theta,s,mem,f) %---------------------------- flowstruct
%FLOWSTRUCT Constructs the structure FL containing the solution.
%  FLOWSTRUCT(THETA,SUBSTANCE,MEMBRANE,FMODEL) constructs a struct FL with the
%  fields FL.INFO, FL.CALC, FL.SOL and FL.FLOW. FLOWSTRUCT sets the fields in
%  FL.INFO and constructs the struct FL.SOL, but does not set any field in
%  FL.SOL. The fields FL.CALC and FL.FLOW remain empty. The struct FL.CALC is
%  constructed in MNUM, FL.FLOW is constructed in FLOW12.
%  FL.info contains
%    .theta
%    .sname
%    .T1
%    .p1
%    .p2
%    .substance
%    .membrane
%    .fmodel
%  FL.sol has the fields
%    .m
%    .T1
%    .p1
%    .q1
%    .T2
%    .T3
%    .Kn2
%    .len
%    .states
%
%  See also FLOW12, FMODEL, MEMBRANE, MNUM, SUBSTANCE.

fl = struct('info',[],'calc',[],'sol',[],'flow',[]);
fl.info = struct('theta',theta,'sname',s.name,'T1',[],'p1',[],'p2',[],...
  'substance',s,'membrane',mem,'fmodel',f,'colors',{{'b','r','g'}});
fl.sol = struct('m',[],'T1',[],'p1',[],'q1',[],'T2',[],'T3',[],'Kn2',[],...
  'len',[],'states',[]);

end %------------------------------------------------------------ end flowstruct
