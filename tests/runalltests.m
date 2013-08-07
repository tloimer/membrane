%RUNALLTESTS A simple reminder to run all tests. Stay in the directory tests/.

[path,name] = fileparts(cd);
if name ~= 'tests', error('Run from directory tests/.'); end
clear('path','name');

if ~exist('substance.m'), addpath('../program'); end

%global VERBOSE;  VERBOSE = 1;
entropy
clear all
frvmCcc
frvmvskap
clear all
testmstack
