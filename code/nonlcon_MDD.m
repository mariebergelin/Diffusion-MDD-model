%% Nonlinear constraints for forward MDD model parameters

function [c,ceq] = nonlcon_MDD(X,Ddata)

% This function set up the nonlinear constraints used for the optimization
% of the diffusion kinetics parameters
%
% The constraint is as follows:
% It is assumed that every diffusion experiment resulted in a
% complete degassing of the isotope. Therefore the modeled total
% fractional release should be 1.

% Constraints are set up such that they satistfy the nonlcon condtitions
% c = 0     % Compute nonlinear inequalities at x
% ceq <= 0  % Compute nonlinear equalities at x
% 
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

% Unpack X
if nargin < 3; plotFlag = 0; end

ndom = (length(X)-2)./2;

% set model parameters
Mtot = X(1);
Ea = X(2);
everything_else = X(3:end);
lnD0aa = everything_else(1:ndom);
D0aa = exp(lnD0aa);
Fx = everything_else(ndom+1:end);

% set model parameters
this_mparam.Mtot = Mtot;
this_mparam.Ea = Ea;
this_mparam.D0aa = D0aa;
this_mparam.Fx = Fx;

% Get observed and modeled data
pred = plot_MDD_forward_model(this_mparam,Ddata,0);

% make sure all Fx values are < 0
%c = [];
c = 0 - min(Fx);

% set constraint such that the modeled total fractional release results in
% a complete degassing

ceq = 1 - pred.Fi(end);

end
