function out = plot_MDD_forward_model(mparam,data,plotFlag,init_Diff)

% plot_MDD_forward_model

% this function plots a forward MDD model for n diffusion domains

% Param has the structure of Ea, D0aa, and Fx (only if more that one domain
% exist). If the sum of FX does not equal 1, an error will appear.

% data should contain the heating schedule, which just needs to contain 
% step of Temp (data.TC) and heating duration (data.tmin). If heating 
% schedule is empty [], then this function will generate a basic increasing
% heating schedule
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

if nargin < 3; plotFlag = 1; end

out = mparam;

ndom = max(size(mparam.D0aa));

% unpack variable parameters
Mtot = mparam.Mtot;
Ea = mparam.Ea;
D0aa = mparam.D0aa;

Fx = mparam.Fx;

% get heating schedule

if isempty(data)
% create a modeled heating schedule
    nstep = 100;
    ti = 30; % min in each heating step
    TK_lim = [50 1500];
    
    TK = 273.15+(TK_lim(1):(TK_lim(2)-TK_lim(1))/(nstep-1):TK_lim(2))';
    t = ones(nstep,1).*ti.*60; % 60 min heating step for every step
    step = (1:size(TK));

else
    % Unpack heating schedule and covert units
    TK = vertcat(data.TC) + 273.15; % convert to degree Kelvin 
    t = vertcat(data.tmin).*60; % convert tmin to tsec
    if isfield(data,'step'); step = data.step;
    else step = (1:size(TK))'; end

end

% also add initial gas loss heating steps in the beginning if applicable
if ~isempty(init_Diff)
    TK = vertcat((init_Diff.TC + 273.15),TK); % convert to degree Kelvin 
    t = vertcat((init_Diff.tmin.*60),t); % convert tmin to tsec
end


% Define constants
R = 0.008314; %gas constant (KJ/(mol K))

% calculate diffusion (D/a^2)
Daa = D0aa.*exp(-Ea./(R.*TK));

% Calculate Cummulative fractional release
% This is taken from Fechtig and Kalbitzer (1966) Equation 4a-c.

% set up equations
Dtaa = Daa.*t;
Dtaa = cumsum(Dtaa); % cummulative sum

Bt = pi.^2.*Dtaa;

Fi = zeros(length(TK),length(D0aa)); %cummulative fractional release
fi = Fi; % fractional release from each heating step

for a = 1:length(D0aa)

    % Cummulative fractional release
        
        Fi(:,a) = 6./(pi^(3/2)).*sqrt(Bt(:,a));         

        % Fechtig and Kalbitzer (1966) Equation 4b
        % Reichenberg (1953) Equation 9
        %i = (min(find(Fi(:,a) > 0.1)):size(Fi(:,a),1));
        ii = find(Fi(:,a) > 0.1);
        i = min(ii):size(Fi(:,a),1);
        Fi(i,a) = (6./(pi^(3/2))).*sqrt(Bt(i,a)) - (3./(pi.^2)).*(Bt(i,a));
       
        % Fechtig and Kalbitzer Equation 4c
        % Reichenberg (1953) Equation 11
        % It should be noted that for long heating steps (e.g Myr) the
        % fractional release resolved from to Eq. 4b/9 decreases after
        % reaching ~0.95, however as noted in Reichenberg (1953), Fi is a
        % mathematical calculable function of Bt (see table 1 in
        % Reichenberg (1953)). at F = 0.99, Bt = 4.11. Any value of Bt >
        % 4.11 results causes the calculated fractional release to decrease
        % and could in some instances reach negative values. To avoid this,
        % find any values of Bt > 4.11 and use Eq. 4c/11.

        %i = (min(find(Fi(:,a) > 0.9 | Bt(:,a) > 1.80)):size(Fi(:,a)))';
        ii = find(Fi(:,a) > 0.9 | Bt(:,a) > 1.80);
        i = min(ii):size(Fi(:,a),1);
        Fi(i,a) = 1 - (6/(pi.^2))*exp(-Bt(i,a));

    Fi(:,a) = Fi(:,a).*Fx(a); % fractional release from each domain
end


for a = 1:length(D0aa)
    % Fractional release from each step
    fi(1,a) = Fi(1,a);
    
    for b = 2:height(Fi)
        fi(b,a) = Fi(b,a)-Fi(b-1,a);
    end
  
end


% Now I want to remove the initial heating steps if they were added
% then scale the release for fi and the calculate Fi

if ~isempty(init_Diff') % this is set up to only handle 3 stages of initial gas loss

    out.Fi_init_Diff = sum(fi(1:3,:),2);

    fi = fi(4:end,:);
    fi = fi./sum(sum(fi,1));
    Fi = cumsum(fi,1);

else out.Fi_init_Diff = zeros(3,1); end

% combined fraction from all domains
Fi_tot = sum(Fi,2); % horizontal summation of cummulative fractional release
fi_tot = sum(fi,2); % horizontal summation each step's fractional release


%%%%%% Here calculate M (predicted measured gas)

M = fi_tot.*Mtot;

% calculate total Daa for multiple domaine model

if ~isempty(init_Diff)
    % now removed added t and TK steps
    t(1:3,:) = [];
    TK(1:3,:) = [];
    Daa(1:3,:) = [];
end

if ndom == 1; Daa_tot = Daa;
else

    % Now calculate the D/a^2 for combined multiple domains

    n = length(step);

    % Fetchig and Kalbitzer Equation 5a
    Daa_tot_a = zeros(n,1);
    
    Fbefore(1,1) = 0; Fbefore(2:n,1) = Fi_tot(1:(end-1),1);
    Fafter = Fi_tot;
    
    Daa_tot_a = (Fafter.^2-Fbefore.^2).*pi./(36.*t);
    
    % Fetchig and Kalbitzer Equation 5b
    Daa_tot_b = zeros(n,1);
    
    Daa_tot_b = (1./((pi.^2).*t)).*( -(pi.*pi./3).*(Fafter-Fbefore) - (2.*pi).*( sqrt(1-(pi./3).*Fafter) - sqrt(1 - (pi./3).*Fbefore) ));
    
    % Fetchig and Kalbitzer Equation 5c
    Daa_tot_c = zeros(n,1);
    
    Daa_tot_c = (1./((pi.^2).*t)).*(log((1-Fbefore)./(1-Fafter)));
    
    % Decide what to use
    
    % sort into apropriate fraction zones
    usea = (Fi_tot <= 0.1 & Fi_tot > 0);
    useb = (Fi_tot > 0.1 & Fi_tot <= 0.9);
    usec = (Fi_tot > 0.9 & Fi_tot <=1);
    
    Daa_tot = usea.*Daa_tot_a + useb.*Daa_tot_b + usec.*Daa_tot_c;
end


% Now do some plotting
if plotFlag == 1

figure(1);

delete(findobj('tag','model'))

    if isempty(findobj('tag','ax_fi')) % if figure does not exist
    
        set(1,'pos',[100   251   1300   300]);

        % set up axes position calculations
        axl = 0.05; axr = 0.1; axsp = 0.05;
        axb = 0.09; axt = 0.05;

        axw = (1-axl-axr-2*axsp)./3;
        axht = (1-axb-axt);

        ax_fi = axes('pos',[axl axb axw axht],'tag','ax_fi'); hold on;
        ax_Fi = axes('pos',[axl+axw+axsp axb axw axht],'tag','ax_Fi'); hold on;
        ax_Arr = axes('pos',[axl+2*axw+2*axsp axb axw axht],'tag','ax_Arr'); hold on;
    end

% Fractional release for each heating step
set(gcf,'currentaxes',findobj('tag','ax_fi'));

    if max(find(D0aa)) > 1
        for a = 1:length(D0aa)
            plot(step,fi(:,a),'color',[1 0 0 0.3],'linewidth',1,'tag','model'); hold on;
        end
    end

    plot(step,fi_tot,'k-','linewidth',1,'DisplayName','Model','tag','model'); hold on;
        
    % add axis labels if not already existing
    if isempty(get(get(gca,'Xlabel'),'string'))
        %legend
        xlabel('Heating step'); ylabel('Fractional Release')
        title('Fractional Release')
    end
    
    
% Cummulative fractional release 
set(gcf,'currentaxes',findobj('tag','ax_Fi'));
  
    if max(find(D0aa)) > 1
        for a = 1:length(D0aa)
            plot(step,Fi(:,a),'color',[1 0 0 0.3],'linewidth',1,'tag','model'); hold on;
        end
    end
    plot(step,Fi_tot,'k','linewidth',1,'DisplayName','Model','tag','model'); hold on;
    

    % add axis labels if not already existing
    if isempty(get(get(gca,'Xlabel'),'string'))
        %legend
        xlabel('Heating step'); ylabel('Cummulative Fractional Release')
        title('Cummulative Fractional Release')
    end
    
   
% Arrhenius plot

set(gcf,'currentaxes',findobj('tag','ax_Arr'));

    for a = 1:length(D0aa)

        mlw = 10*Fx(a);
        if mlw < 0.5; mlw = 0.5; end
        % plot all domains line
        i = find(TK == min(TK));
        ii = find(TK == max(TK));
        h = plot([1e4./TK(i(1)) 1e4./TK(ii(end))],[log(Daa(i(1),a)) log(Daa(ii(end),a))],'--','color',[1 0 0 0.3],'linewidth',mlw,'DisplayName',['Model Domains' sprintf('%0.0f',a)],'tag','model');
        uistack(h,'bottom')
    end

    plot((1e4./TK),log(Daa_tot),'k.','markersize',12,'tag','model');
    
    % add axis labels if not already existing
    if isempty(get(get(gca,'Xlabel'),'string'))
        legend
        xlabel('10^4/T [K^-^1]'); ylabel('ln(D/a^2) [ln(s^-^1)]')
        title('Arrhenius Plot')
    end
end

out.Fx = Fx;
out.n = max(step);
out.step = step;
out.TC = data.TC;
out.tmin = data.tmin;
out.fidom = fi;
out.fi = fi_tot;

if D0aa > 1
out.Fidom = Fi; end

out.Fi = Fi_tot;
out.Daadom = Daa;
out.Daa = Daa_tot;

out.M = M;

end