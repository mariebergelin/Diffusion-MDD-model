function out = get_MDD_misfit(mparam,Ddata,plotFlag,init_Diff)

% Various parameters to match observed and predicted with.
% .M = gas released from each step
% .fi = fractional release from each step
% .Fi = cummulative fractional release
% .Daa = Diffusion

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

% if plotFlag == 1; delete(findobj('tag','model')); end

% get data from Ddata diffusion experiment
obs = plot_Ddata(Ddata,0);

% get predicted release data from forward model
pred = plot_MDD_forward_model(mparam,Ddata,plotFlag,init_Diff);

% Get which steps to include in the misfit calculations
mask = find(obs.use > 0);

% degrees of Fredom
nparam = 2+(length(mparam.D0aa).*2); % number of parameters
df = length(mask)-nparam; % degree of fredom

% Now compute misfits

% Percent Fraction misfit from Gorin et al. (2024)
Mi = abs(obs.fi(mask)-pred.fi(mask))./obs.fi(mask);
Mi(find(isinf(Mi))) = 0; Mi(find(isnan(Mi))) = 1; % find Inf numbers and make zero as this will provide problems
M = sum(Mi);

% % Other statistical misfit options to mess around with

% % Chi^2
% Mi = abs(obs.M(mask)-pred.M(mask))./obs.delM(mask);
% M = sum(Mi.^2);

% % Tremblay et al. (2017)
% 
% Mi = abs(obs.Fi(mask)-pred.Fi(mask));
% M = sum(Mi);

% % absolute residual
% Mi = (obs.fi(mask)-pred.fi(mask)).^2;
% M = sum(Mi);

% penalize complete degassing of sample before the last step
%      i = find(pred.Fi(mask) > 0.9999);
%      Mi(i(1:end-2)) = 0.1;
%     Mi(find(Mi > 0.98 & Mi <= 1)) = 1000;

% % misfit from Daa
% Mi = abs((obs.Daa(mask) - pred.Daa(mask))./obs.Daa(mask));
% Mi(find(Mi > 0.98 & Mi <= 1)) = 100;
% M = sum((Mi(~(isnan(Mi) | isinf(Mi)))).^2);


% Mi = abs(obs.M(mask)-pred.M(mask));
% Mi = abs(obs.fi(mask)-pred.fi(mask))./obs.fi(mask);

% M = sum(Mi.^2);
% M = sum(Mi);


% Make a plot of the misfit in each step

% if plotFlag == 1
%     % plot misfit
%     
% figure(2);
% 
%     if isempty(get(get(gca,'Xlabel'),'string')) % if figure does not exist
%         clf;
%         % set-up new figure window
%         set(2,'pos',[300   15   800   350]);
%     
%         % set up axes position calculations
%         axl = 0.08; axr = 0.05; axsp = 0.05;
%         axb = 0.09; axt = 0.05;
%         
%         axw = (1-axl-axr-axsp)./2;
%         axht = (1-axb-axt);
%         
%         ax_MS1 = axes('pos',[axl axb axw axht],'tag','ax_Mi'); hold on;
%         ax_MS2 = axes('pos',[axl+axw+axsp axb axw axht],'tag','ax_M'); hold on;
% 
%     end
%     
% set(gcf,'currentaxes',findobj('tag','ax_Mi'));
%     
%     set(findobj('tag','Mi'),'color',[1 0.85 0.85]); % Change the color of previous plots
%     plot(obs.step(mask),Mi,'-r','tag','Mi'); hold on;
%     
%     set(gca,'ylim',[0.01 10e5],'yscale','log')
%     set(gca,'xlim',[0 max(obs.step)])
%     xlabel('step')
%     ylabel('Misfit')
% 
% set(gcf,'currentaxes',findobj('tag','ax_M'));
%     MF_count = get(2,'UserData');
% 
%     %if isempty(findobj('tag','misfit')); MF_count = 0; end
%     if isempty(MF_count); MF_count = 0; end
%     
%     if MF_count > 1
%         ax = get(findobj('tag','misfit'),'color');
%         if all(ax{1} == [1 0.8500 0.8500]); MF_count = 0; end
%             % This changes the previous MC iteration plotted in red data to faded red
%     end
%     
%     MF_count = MF_count + 1;
%     set(2,'UserData', MF_count);
%     
%     %set(gca,'xlim',[0 MF_count+1],'ylim',[0.1 1000],'yscale','log')
%     set(gca,'ylim',[0.01 10e6],'yscale','log')
%     set(gca,'xlim',[0 300])
% 
%     % plot every M
%     plot(MF_count,M,'.r','markersize',10,'tag','misfit'); hold on
%     xlabel('iteration')



% return misfit value
out = M; end
