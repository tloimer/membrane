% Flow through porous membranes
% Version asym 02-Apr-2015
%
% Setup.
%   substance     - Material properties of a substance.
%   fmodel        - Two-phase flow model.
%   membrane      - A homogeneous membrane.
%   mstackstruct  - A stack of membranes, each consisting of several layers.
%   readdata      - Read data from a file.
%
% Compute.
%   mnumadiabat   - Mass flux for adiabatic flow through layered membranes.
%   asym          - Upstream pressure for adiabatic flow through a MSTACKSTRUCT.
%   mgaseous      - Mass flux for isothermal gaseous flow.
%   misotherm     - Isothermal mass flux through a homogeneous membrane.
%   mlinear       - Mass flux from linear theory through a homogeneous membrane.
%   calctaubeta   - Calculate tau and beta.
%   See also functions in MSTACKSTRUCT, see below.
%
% Plot.
%   beginpgfplot  - Write the begin of a pgfplot.
%   addcoords     - Add data to a pgfplot.
%   endpgfplot    - Close a pgfplot.
%   pTplot        - Plot p-T diagram.
%   pTzplots      - Plot p-z, T-z and, optionally, q-z diagrams.
%   Tsplot        - Plot a T-s Diagram.
%   See also functions in MSTACKSTRUCT, see below.
%
% Utilities.
%   downstreamstate - Return a STATE struct.
%   findinterval    - Find an interval in which a function changes sign.
%   findzero        - Find the zero of a function.
%   flowsetup       - Setup flow properties.
%   solverstruct    - Construct a solver struct.
%
% Data.
%   mem_data - Properties of the Hermsdorf membranes.
%
% See also mstackstruct>mfluxliquid, mstackstruct>mfluxviscous,
%          mstackstruct>mfluxknudsen, mstackstruct>printsetup,
%          mstackstruct>plotsolution, mstackstruct>printsolution.
