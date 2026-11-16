function out = objective_MDD_model(X,Ddata,plotFlag,init_Diff)

% This function is a wrapper for the get_MDD_misfit that allow it to be
% used with the MATLAB optimizer fminsearch()
%
% This code is provided for the sole purpose of evaluation of an
% accompanying scientific paper preprint entitled ’Diffusion kinetics of
% 3He in pyroxene and plagioclase and applications to cosmogenic exposure
% dating and paleothermometry in mafic rocks’ by Marie Bergelin and three
% other authors. It is not licensed for any other reuse or redistribution.
%
% Written by Marie Bergelin 
% Contact: mbergelin@bgc.org
% Last modified: 2025.02.27

if nargin < 4; init_Diff = []; end

if nargin < 3; plotFlag = 0; end

ndom = (length(X)-2)./2;

% set model parameters
this_mparam.Mtot = X(1);
this_mparam.Ea = X(2);
everything_else = X(3:end);

lnD0aa = everything_else(1:ndom);
this_mparam.D0aa = exp(lnD0aa);

this_mparam.Fx = (everything_else(ndom+1:end));

out = get_MDD_misfit(this_mparam,Ddata,plotFlag,init_Diff); end