% Generate a subplot from MDDdata results

% This script generate Fig. 6 and Fig. 7 in Bergelin et al., "Diffusion
% kinetics of 3He in pyroxene and plagioclase and applications to
% cosmogenic exposure dating and paleothermometry in mafic rocks"

% For plotting, this script requires the filename for the diffusion
% experiment text file and the filename for the MDDdata results matlab structure

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

% list Ddata file name in order for plotting
fndata = {'LABCO-Px-Db.txt','LABCO-Px-Dc.txt','ROB-Px-Da.txt','GEM-CPx-Da.txt','GEM-OPx-Da.txt'};
fnMDD = {'LABCO-Px-Db_final.mat','LABCO-Px-Dc_final.mat','ROB-Px-Da_final.mat','GEM-CPx-Da_final.mat','GEM-OPx-Da_final.mat'};

% Unpack Ddata structures
for a = 1:length(fndata)
    
    thisDdata = unpack_Ddata(char(fndata(a)));
 
    Ddata{a} = thisDdata;
    Ddata{a} = plot_Ddata(Ddata{a},0);

    thisMDDdata = importdata(char(fnMDD(a)));
    MDDdata{a} = thisMDDdata;

end

% get color order and reorder
% i = [1,2,3,4,5,7,9,10];
% col = colororder("gem12"); col = col(i,:);

greens = [0.2745    0.4196    0.0824
    0.4235    0.6392    0.1412
    0.5922    0.7608    0.3725];

purples = [0.4941    0.1843    0.5569
    0.6745    0.5176    0.7686];

col = [greens; purples];

% Make faded color
for a = 1:height(col)
    thiscol = col(a,:);
    % make color faded
    [thiscols, i] = sort(thiscol);
    thiscol0 = thiscol + (1-thiscol)./1.3;
    thiscol0(i(3)) = 0.95;
    thiscolr = thiscol(i(2))/thiscol(i(1));
    thiscol0(i(1)) = thiscol0(i(2))/1.2; %(col_r/3);
    col0(a,:) = thiscol0;
end

fig_ann = {'a)','b)','c)','a)','b)'};

for a = 1:length(fnMDD)
    
    thisDdata = Ddata{a};
    thisMDDdata = MDDdata{a};

    msz = 8; % marker size
    lw = 1; % line width


    Tplot = 1e4./(thisDdata.TC + 273.15); % x values for plotting arrhenius

    
    figure(a);
    
    set(gcf,'pos',[100   450   480   300]);
    annotation('textbox',[0.01 0.92 0.1 0.1],'String',fig_ann(a),'FontSize',16,'FitBoxToText','on','EdgeColor','none');

    
    % set up axes position calculations
        axl = 0.1; axr = 0.02; 
        axvsp = 0.15; axhsp = 0.11;
        axb = 0.14; axt = 0.09;

        axw1 = (1-axl-axr-axhsp)./1.5;
        axht1 = (1-axb-axt);
        axw2 = (1-axl-axr-axhsp-axw1);
        axht2 = (1-axb-axt-axvsp)./2;

        axoff = 0.1; axb0 = 0.1;

    % Arrhenius plot
    ax_Arr = axes('pos',[axl axb+axoff axw1 axht1-axoff],'box','on','FontSize',12); hold on;
       
        plot(Tplot,log(thisDdata.Daa),'o','color',col0(a,:),'markersize',msz,'linewidth',lw);
    
        % Assemble parameters for fitting
        fitOK = find(thisDdata.use > 0);
        fitDaa = thisDdata.Daa(fitOK);
        fitdelDaa = thisDdata.delDaa(fitOK);
        fitT = Tplot(fitOK);
    
        % Ln error propagation and plotting
        for b = 1:length(fitOK)
            fitdellnDaa(b) = fitdelDaa(b)./fitDaa(b);
            % yy1 = fitDaa(a) - fitdelDaa(a);
            %     if yy1 < 0; yy1 = exp(-35); end
    
            xx = [fitT(b) fitT(b)];
            %yy = [log(yy1) log(fitDaa(a)+fitdelDaa(a))];
            yy = [(log(fitDaa(b))-fitdellnDaa(b)) (log(fitDaa(b))+fitdellnDaa(b))];
            plot(xx,yy,'-','color',col(a,:),'linewidth',lw);
        end
    
        % plot data listed for use
        %plot(fitT,log(fitDaa),'.','MarkerSize',25,'color',col_0,'tag','data'); hold on;
        plot(fitT,log(fitDaa),'o','MarkerSize',msz,'linewidth',lw,'color',col(a,:),'markerfacecolor',col0(a,:),'tag','data'); hold on;
    
    
        % plot model results
        for b = 1:length(thisMDDdata.D0aa)
    
            % plot all domains line
            %h = plot([Tplot(1) Tplot(end)],[log(MDD_results.Daadom(1,a)) log(MDD_results.Daadom(end,a))],'--','color',[1 0 0 0.3],'linewidth',10*MDDdata.Fx(a),'DisplayName',['Model Domains' sprintf('%0.0f',a)]);
            %line with
            mlw = 10*thisMDDdata.Fx(b);
            if mlw < 0.6; mlw = 0.6; end
            i = find(Tplot == max(Tplot));
            ii = find(Tplot == min(Tplot));
            h = plot([Tplot(i(1)) Tplot(ii(end))],[log(thisMDDdata.Daadom(i(1),b)) log(thisMDDdata.Daadom(ii(end),b))],'--','color',[0 0 0 0.25],'linewidth',mlw,'DisplayName',['Model Domains' sprintf('%0.0f',a)]);
            uistack(h,'bottom')
        end
    
        plot(Tplot,log(thisMDDdata.Daa),'k.','markersize',12,'tag','model');
    
        title(thisDdata.samplename)
        set(gca,'xlim',[5 35],'ylim',[-40 -5]);
    
        xlabel('10^4/T(K)');
        ylabel('ln(D/a^2)');
    

        ax0 = axes('pos',[axl axb0 axw1 axoff/4],'box','off','color','none');

        tickts = fliplr([25 50 100 200 300 400 500 600 700 800 900 1000 1100 1200]+273.15);
        tickxs = 10000./tickts;
        tickls = {'1200','','','','','','600','','400','300','200','100','50','25'};
        set(gca,'xlim',get(ax_Arr,'xlim'),'xtick',tickxs,'xticklabel',tickls,'ytick',[],'Yticklabel',[]);
        xtickangle(0);
        xlabel('T (°C)');
        set(gca,'fontsize',10);
    
    % Fractional release plot
    ax_Fi = axes('pos',[axl+axw1+axhsp axb+axht2+axvsp axw2 axht2],'box','on','FontSize',10); hold on;
    
        plot(thisDdata.step,thisDdata.fi,'-','color',col0(a,:)); hold on;
    
        plot(thisDdata.step,thisDdata.fi,'o','MarkerSize',msz,'linewidth',lw,'color',col(a,:));
        plot(thisDdata.step(fitOK),thisDdata.fi(fitOK),'o','MarkerSize',msz,'linewidth',lw,'color',col(a,:),'markerfacecolor',col0(a,:));
    
        % Plot model
    
        plot(thisMDDdata.step,thisMDDdata.fi,'-k');
        plot(thisMDDdata.step,thisMDDdata.fi,'.k','MarkerSize',10);
    
        xlabel('heating step');
        ylabel('fractional gas release');
    
    % residual plot
    ax_res = axes('pos',[axl+axw1+axhsp axb axw2 axht2],'box','on','FontSize',12);
    
        % get Ddata residual from domain number 1
        resDdata = log(thisDdata.Daa)- (thisMDDdata.lnD0aa(1)-(thisMDDdata.Ea/83.14)*(1e4./(thisDdata.TC+273.15)));
        % get MDDDdata residual from domain number 1
        resMDDdata = log(thisMDDdata.Daa)- (thisMDDdata.lnD0aa(1)-(thisMDDdata.Ea/83.14)*(1e4./(thisDdata.TC+273.15)));
        
        % calculate residual to be plotted
        res = (log(thisDdata.Daa) - log(thisMDDdata.Daa))./log(thisDdata.delDaa);
        
        % cummulative fractional release
        cumFi = thisDdata.Fi;
        
        plot(cumFi,res,'o','MarkerSize',8,'linewidth',1,'color',col(a,:),'markerfacecolor',col0(a,:)); hold on;
        
        %plot(cumFi,res,'.','color',[0 0 0]);
        
        % plot(cumFi,resDdata,'o','MarkerSize',8,'linewidth',1,'color',col(a,:),'markerfacecolor',col0(a,:)); hold on;
        % plot(cumFi,resDdata,'--','color',col0(a,:));
        % 
        % plot(cumFi,resMDDdata,'.k','markersize',12);
        % plot(cumFi,resMDDdata,'--k');
    
        xlabel('Cummulative gas release');
        ylabel('residual (ln(s^-^2))');

        xlim([-0.02 1.02]);
        ylim(get(gca,'Ylim')+ 0.02*[-1 1]);

        plot(get(gca,'Xlim'),[0 0],'--','color',[0.3 0.3 0.3],'linewidth',0.5); hold on;



end
