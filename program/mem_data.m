function [mf,mb] = mem_data(name)
%MEM_DATA   Properties of the Hermsdorf membranes.
%  MEM_DATA(NAME) returns the memstackstruct with the properties of a
%  Hermsdorf membrane. NAME can be 'hermsdorf2001' or 'hermsdorf2014'.
%
%  [MF,MB] = MEM_DATA(NAME) returns the memstackstructs MF and MB for
%  forward and backward flow.
%
%  If the global variable VERBOSE is greater than zero, the membrane
%  properties are displayed.
%
%  See also MSTACKSTRUCT and FMODEL.

global VERBOSE

f = fmodel('emt');	% Thermal conductivity model
theta = 0;		% Contact angle

switch (name)
case 'hermsdorf2014'
    % From the support to the separation layer
    %		 pore diameter	epsilon	km		tau	beta	L
    l1 = membrane(2.5e-6,	0.6,	30,	'tube',	4.8,	8.1,	1e-3);
    l2 = membrane(0.8e-6,	0.5,	30,	'tube',	3.05,	8.1,	25e-6);
    l3 = membrane(0.4e-6,	0.5,	30,	'tube',	2.66,	8.1,	25e-6);
    l4 = membrane(0.1e-6,	0.4,	30,	'tube',	1.6,	8.1,	25e-6);
    l5 = membrane(30e-9,	0.3,	11.8,	'tube',	3,	8.1,	5e-6);

case 'hermsdorf2001'
    l1 = membrane(3.07e-6,	0.4,	1,	'tube',	3.2,	9.0541,	1.5e-3);
							% eps/tau = 0.125
    l2 = membrane(1.84e-6,	0.5,	1,	'tube',	3.05,	9.0541,	25e-6);
							% eps/tau = 0.164
    l3 = membrane( 191e-9,	0.5,	1,	'tube',	2.66,	9.0541,	25e-6);
							% eps/tau = 0.188
    l4 = membrane(  76e-9,	0.6,	1,	'tube',	1.51,	9.0541,	25e-6);
							% eps/tau = 0.397
    l5 = membrane ( 16e-9,	0.8,	11.8,	'tube',	1.08,	8.1,	300e-9);
							% eps/tau = 0.741
otherwise
    error('Use ''hermsdorf2001'' or ''hermsdorf2014'' as argument to %s.',...
	mfilename);
end

mf = mstackstruct(theta, {{l5 l4 l3 l2 l1}},f);
if nargout > 1
    mb = mstackstruct(theta, {{l1 l2 l3 l4 l5}},f);
end

if VERBOSE > 0
    mf.printsetup(mf);
end
