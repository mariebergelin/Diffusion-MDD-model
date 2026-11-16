function out = plot_Ddata(data,plotFlag,col)

% This function calculates the fractional realease and diffusion for each
% step in a diffusion experiment and plots the fractional release,
% cumulative fractional release and the Arrhenius plot, and returns all the
% calculated parameters in a structure file for the diffusion experiment.
%
% The input for this function is the data structure obtained from
% unpack_Ddata. If plotting multiple experiments this function can be used
% in a loop which provide different colors for each data experiment (see
% example below). If no input is included, this function will ask for a
% single data .txt files to be selected as input.
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
%
%
%%%%%%%% Example for use in a script:
%
% % Mulitple files can be plotted
% [fn,pn] = uigetfile('*.txt','All Files (*.*)','MultiSelect','on');
%
% plotFlag = 1;   % Plot data
% 
% % Plot diffusion data
% if isa(fn,'char')
%     Ddata = unpack_Ddata(fn);
%     Ddata = plot_Ddata(Ddata,dplotFlag);
% else
%     colororder = get(gca,'ColorOrder');
%     
%     for a = 1:length(fn)
%     
%         col = colororder(a,:);
%         thisdata = unpack_Ddata(char(fn(a)));
%      
%         Ddata{a} = thisdata;
%         Ddata{a} = plot_Ddata(Ddata{a},plotFlag,col);
% 
%     end
% end
%
%%%%%%%%%%%%%


if nargin < 3 
    col_1 = [0 0.45 0.74]; % blue
    col_0 = [0.73 0.83 0.89]; % faded blue (data not used)
elseif exist('col','var')
    col_1 = col;

    % Make the input color faded (data not used)
    [col1s, i] = sort(col_1);
    col_0 = col_1 + (1-col_1)./1.3;
    col_0(i(3)) = 0.95;
    col_r = col_1(i(2))/col_1(i(1));
    col_0(i(1)) = col_0(i(2))/1.2; %(col_r/3);

end 

if nargin < 2; plotFlag = 0; end

if nargin < 1
    plotFlag = 1;
    [fn,pn] = uigetfile('*.txt');
    data = unpack_Ddata(fn);
end

% unpack data
TC = data.TC;       % heating temperature in degree C
tmin = data.tmin;   % heating duration in minutes
M = data.M;         % gas amount released, unite doesn't matter
delM = data.delM;   % uncertainty in gas amount released 
n = length(M);

if isfield(data, 'use')
else data.use = ones(data.n,1); end

% Rationalize everything to column vectors
if size(data.TC,1) < size(data.TC,2)
    data.TC = data.TC';
end
if size(data.tmin,1) < size(data.tmin,2)
    data.tmin = data.tmin';
end
if size(data.M,1) < size(data.M,2)
    data.M = data.M';
end
if size(data.delM,1) < size(data.delM,2)
    data.delM = data.delM';
end

use_mask = find(data.use > 0); 

% convert units
TK = 273.15+TC; 
tsec = tmin.*60;
cumtsec = cumsum(tsec);
Tplot = 1e4./TK;

% sums
Si = cumsum(M);
delSi = cumsum(delM);
S = max(Si);
delS = sqrt(sum(delM.^2));

% fraction released during degassing

fi = M./S; % fraction released during each step

% Uncertainty calculated from Ginster and Reiners (2018) Equation 17 
for a = 1:n
    delfi(a,1) = (1/S).* sqrt(((1-(M(a)./S)).*delM(a)).^2+(M(a)./S).^2.*(sum(delM(1:a).^2)+(sum(delM(a:end).^2))));
end

Fi = Si./S; % cummulative fractional release

% Uncertainty calculated from Ginster and Reiners (2018) Equation 17 
delFi = zeros(n,1);

for a = 1:n-1
    delFi(a,1) = (1/S).* sqrt((1-Fi(a)).^2.*sum(delM(1:a).^2)+(Fi(a)).^2.*sum(delM(a+1:end).^2));
end

% Now calculate the Da^2

% Ginster and Reiners (2018) Table 4, Equation 27 - 34
Daa = zeros(n,1); delDaa = Daa; uncert_b = delDaa;

% The first step

if Fi(1) <= 0.85

    Daa(1) = 1./(pi.^2.*tsec(1)).* (2.*pi-(pi.^2./3)*Fi(1)-2.*pi.*sqrt(1-(pi./3)*Fi(1)));
    delDaa(1) = (1./(3*tsec(1).*S)).*((1./(sqrt(1-(pi./3*Fi(1)))))-1)*sqrt(((1-Fi(1)).*delM(1)).^2+Fi(1).^2*sum(delM(2:end).^2));
else 
    % if Fi(1) < 0.85
    Daa(1) = -1./(pi.^2.*tsec(1)).*log((1-Fi(1)).*pi.^2./6);         
    delDaa(1) = 1./(pi.^2*tsec(1).*S).*sqrt(delM(1).^2+(Fi(1)./(1-Fi(1))).^2.*sum(delM(2:end).^2));
end

% all other steps
for a = 2:n
    
    % dummy eqations used in repetitive uncertainty calculations
    AFi = 1/sqrt(1-(pi./3*Fi(a))); BFi = 1/sqrt(1-(pi./3*Fi(a-1)));

    if Fi(a) <= 0.85

        Daa(a) =  1./(pi.^2*tsec(a)).*(((-pi.^2./3)*(Fi(a)-Fi(a-1)))-(2.*pi.*(sqrt(1-(pi./3).*Fi(a))-sqrt(1-(pi./3).*Fi(a-1)))));
        
        delDaa(a) = 1./(3*tsec(a).*S).* sqrt(((1-Fi(a)).* (AFi-1) + (1-Fi(a-1)).* (1-BFi)).^2.* sum(delM(1:a-1).^2) ...
                                                + ((1-Fi(a)).*(AFi-1) - Fi(a-1).*(1-BFi)).^2   .*  delM(a).^2  ...
                                                + (Fi(a).*(1-AFi) - Fi(a-1)*(1-BFi)).^2   .*  sum(delM(a+1:end).^2));
    else
        Daa(a) =  -1./(pi.^2*tsec(a)).*log((1-Fi(a))/(1-Fi(a-1)));
        
        delDaa(a) = 1./(pi.^2*tsec(a)*S).* sqrt((1+Fi(a-1)./(1-Fi(a-1))).^2 .* delM(a).^2   + ((Fi(a-1)./(1-Fi(a-1)))-(Fi(a)./(1-Fi(a)))).^2.* sum(delM(a+1:end).^2)                   ) ;
        uncert_b(a) = (1/(pi.^2*(cumtsec(a)-cumtsec(a-1))*sum(M))) * sqrt((1+(Fi(a-1)/(1-Fi(a-1)))).^2 * delM(a).^2 + (Fi(a-1)/(1-Fi(a-1)) - Fi(a)/(1-Fi(a))).^2 * sum(delM(a+1:end).^2));
    end
end

%% Now make some plots

if plotFlag == 1

    sn = strrep(data.samplename, '_', ' ');

figure(1);
    
    if isempty(get(get(gca,'Xlabel'),'string')) % if figure does not exist
        clf;
        f1 = figure(1);
    
        set(gcf,'pos',[100   450   1300   300])

        % set up axes position calculations
        axl = 0.05; axr = 0.1; axsp = 0.05;
        axb = 0.09; axt = 0.05;
%        axb = 0.17; axt = 0.07;

        axw = (1-axl-axr-2*axsp)./3;
        axht = (1-axb-axt);

        ax_fi = axes('pos',[axl axb axw axht],'tag','ax_fi'); hold on;
        ax_Fi = axes('pos',[axl+axw+axsp axb axw axht],'tag','ax_Fi'); hold on;
        ax_Arr = axes('pos',[axl+2*axw+2*axsp axb axw axht],'tag','ax_Arr'); hold on;
    else axr = 0.1; end
    

% Fractional release for each heating step
set(gcf,'currentaxes',findobj('tag','ax_fi'));

    % plot error shade
    %fill([data.step',fliplr(data.step')],[(fi+delfi)',fliplr((fi-delfi)')],col_0,'edgecolor','none');

    plot(data.step,fi,'-','color',col_0,'linewidth',1,'HandleVisibility','off','tag','data'); hold on;
    plot(data.step(use_mask),fi(use_mask),'-','color',col_1,'linewidth',2,'DisplayName',sn,'tag','data'); hold on;
    
    xlabel('Heating step'); ylabel('Fractional Release')
    title('Fractional Release')
    l = legend;
    set(l,'pos',[1.01-axr 0.8 0.08 0.08]);
    max(get(gca,'ylim'));
    set(gca,'ylim', [0 0.3],'xlim',[1 50]);
    %set(gca,'ylim', [0 0.2],'xlim',[1 50]);

    if max(get(gca,'xlim')) < n
        set(gca,'xlim',[1 n]);
    end
    box on


% Cummulative fractional release 
set(gcf,'currentaxes',findobj('tag','ax_Fi'));

    % plot error shade
    %fill([data.step',fliplr(data.step')],[(Fi+delFi)',fliplr((Fi-delFi)')],col_0,'edgecolor','none');

    plot(data.step,Fi,'-','color',col_0,'linewidth',1,'tag','data'); hold on;
    plot(data.step(use_mask),Fi(use_mask),'-','color',col_1,'linewidth',2,'tag','data'); hold on;
    
    xlabel('Heating step'); ylabel('Cummulative Fractional Release')
    title('Cummulative Fractional Release')
    
    set(gca,'xlim',[1 50],'ylim',[0 1]);

    if max(get(gca,'xlim')) < n
        set(gca,'xlim',[1 n],'ylim',[0 1]);
    end
    box on

% Arrhenius plot
 set(gcf,'currentaxes',findobj('tag','ax_Arr'));
    msz = 8;
    lw = 1;
    box on
    
    % plot all Ddata
    plot(Tplot,log(Daa),'o','color',col_0,'markersize',msz,'linewidth',lw,'tag','data'); hold on;    
   
    % Assemble parameters for fitting
    fitOK = find(data.use > 0);
    fitDaa = Daa(fitOK);
    fitdelDaa = delDaa(fitOK);
    fitT = Tplot(fitOK);
    
    % Ln error propagation and plotting
    for a = 1:length(fitOK)
        fitdellnDaa(a) = fitdelDaa(a)./fitDaa(a);

        xx = [fitT(a) fitT(a)];
        yy = [(log(fitDaa(a))-fitdellnDaa(a)) (log(fitDaa(a))+fitdellnDaa(a))];
        plot(xx,yy,'-','color',col_1,'linewidth',lw);
    end

    % plot data listed for use
     plot(fitT,log(fitDaa),'o','MarkerSize',msz,'linewidth',lw,'color',col_1,'markerfacecolor',col_0,'tag','data'); hold on;

    xlabel('10^4/T [K^-^1]'); ylabel('ln(D/a^2) [ln(s^-^1)]')
    title('Arrhenius Plot')
    set(gca,'xlim',[5 35],'ylim',[-40 0]);

end

out = data;
out.fi = fi;
out.delfi = delfi;
out.Fi = Fi;
out.Daa = Daa;
out.delDaa = delDaa;

end
