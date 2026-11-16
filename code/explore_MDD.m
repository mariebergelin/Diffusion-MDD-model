% This script loads a single diffusion experiment data file and calculate 
% the Do and Ea for a single domain

% This script compares diffusion experiment results with a
% given set of diffusion kinetics either defined within the code or
% selected as a m-file. Further, this script allows you to explore the
% affect that the initial irradiation and/or storage loss conditions has on
% the apparent diffusivity.

% When running: 
% 1. a file selection window will pop-up, here you will be able to select
% the text file containing diffusion experiment results that you would like
% to fit the MDD model to. If multiple txt-files are selected for plotting,
% the forward MDD model is calculated only for the first selected data.
% 2. If diffusion kinetics parameters are defined within the script, then
% the 'InitDiff_GUI' will appear. This gui allows you to input the time and
% temperature of initil gas loss conditions to be include in the forward
% model. selceting 'cancel' or inputting 'zero' values will exclude initial
% gas loss from the forward MDD model. 
% 3. If no diffusion kinetic parameters are defined (hence, commented out),
% then a window will pop-up to selected m-files containing MDD diffusion
% kinetic results. This m-file is the MDDdata optimization results saved
% when running the 'optimizer_MDD_Domain.m' script.
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


clear all; close all;

% Mulitple files can be
[fn,pn] = uigetfile('*.txt','All Files (*.*)','MultiSelect','on');

dplotFlag = 1;   % Plot Ddata
mplotFlag = 1;   % Plot Model

% Plot diffusion data
if isa(fn,'char')
    Ddata = unpack_Ddata(fn);
    Ddata = plot_Ddata(Ddata,dplotFlag);
else
    colororder = get(gca,'ColorOrder');
    
    for a = 1:length(fn)
    
        col = colororder(a,:);
        thisdata = unpack_Ddata(char(fn(a)));
     
        Ddata{a} = thisdata;
        Ddata{a} = plot_Ddata(Ddata{a},dplotFlag,col);

    end
end

%% Set diffusion model parameters

clear Ea

% % Define diffusion kinetic parameters
% % uncomment to explore different MDD diffusion kinetics
% Ea = 80; % Activation Energy (KJ/mol)
% lnD0aa = [8 6 3]; % Diffusion (s^-1) for domain [1 2 3 etc.]
% Fx = [0.70 0.27 0.3]; % fraction release from domain [1 2 etc. ]
% 
% % if multiple diffusion experiment are selected plot the forward model for
% % Ddata{1}
% if length(Ddata) == 1; Mtot = sum(Ddata.M); else Mtot = sum(Ddata{1}.M); end


% or get input from MDDdata results
if ~exist('Ea','var')

    [fn pn] = uigetfile('*.mat');
    MDDdata = importdata(fn);

    Ea = MDDdata.Ea;
    lnD0aa = MDDdata.lnD0aa;
    Fx = MDDdata.Fx;
    Mtot = MDDdata.Mtot;
    init_Diff = MDDdata.init_Diff;

else
    init_Diff = initDiff_GUI(); % Ask if user want to include initial diffusion steps
end

%pack diffusion kinetics into a structure for plotting
mparam.Mtot = Mtot;
mparam.Ea = Ea;
mparam.D0aa = exp(lnD0aa);
mparam.Fx = Fx;

% Plot Forward model
if iscell(Ddata)
    MDD_results = plot_MDD_forward_model(mparam,Ddata{1,1},mplotFlag,init_Diff);
else
    MDD_results = plot_MDD_forward_model(mparam,Ddata,mplotFlag,init_Diff);
end
