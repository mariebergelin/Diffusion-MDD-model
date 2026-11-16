% Forward MDD Optimization model

% This scripts loads a diffusion experiment text file and then use a 
% forward diffusion model to optimize for the diffusion kinetetic 
% parameters: the Activation Energy (Ea, KJ/mol), Diffusion (lnDaa) and the
% fractional release (Fx) through a global search using the fmincon function 
% 
% This optimization script accepts a diffusion experiment Ddata .txt file
% that consist of data lines with 5-6 columns. The first three cols are
% step number, Temp (Deg C), and time (min). Subsequent pairs of columns
% are amount of diffusant and uncertainty therein, ideally the unit should
% be in moles. The last column (col #6) is the data use for modelfitting
% and is optional. If exist, this column should have zeroes and ones,
% indicate which data step to be used in the analysis. If absent, the
% script will automatically generate a generic data use column of ones, which
% includes all steps for the analysis.
% 
% Before running, make the appropriate changes in the script. 
% 1. Make sure correct data is loaded, or window to select data pops up. 
% 2. Set initial starting values and boundaries for the model. These can
% either be the best-guess values or random selected values within the
% accepted boundary. Unless specified the default is: Ea = 65, lnD0aa is
% randomly seleceted to decend with domain number from 17 and 3, and FX is
% randomly selected to decend with domain number between 1 and 0.
% 3. Set plotFlag (1 = plot data and model fitting, 0 = suppress plotting)
% turning OFF plotflag will make the model much much much faster, and
% recommended for all optimization runs. turning ON plotflag is recommened
% for visual and debugging.
% 4. set the number of total optimization interations (numits) for each
% number of domains optimized for. Default numits = 20;

% To run this script, the following functions are needed:
% unpack_Ddata.m
% initDiff_GUI.m
% plot_Ddata.m
% objective_MDD_model.m
% get_MDD_misfit.m
% plot_MDD_forward_model.m
% nonlcon_MDD.m

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

% Get diffusion data
[fn, pn] = uigetfile('*.txt');

% define diffusion prior to experiment
init_Diff = initDiff_GUI;

%%

% unpack the data txt file
Ddata = unpack_Ddata(fn);

% starting parameter for Mtot
Mtot = sum(Ddata.M); % atoms
delMtot = sqrt(sum(Ddata.delM.^2));

% Upper/Lower Boundaries
Mtotlim = [(Mtot-2*delMtot) (Mtot+2*delMtot)];
Ealim = [0 300];
lnD0aalim = [-10 55];
Fxlim = [1e-9 1];

gs = GlobalSearch('Display','off','NumTrialPoints',100000);
%gs = GlobalSearch('Display','off','NumTrialPoints',100000,'DistanceThresholdFactor',0.2,'PenaltyThresholdFactor',0.8, 'XTolerance',1e-9,'BasinRadiusFactor',0.8,'MaxWaitCycle',40);
ms = MultiStart(gs,'UseParallel',true,'Display','off','FunctionTolerance',1e-12,'XTolerance',1e-12);

% number of multistart
nms = 100; % 100*length(X0);

% minimum and maximum number of domains to optimize for
minndom = 2;
maxndom = 10;

% initiat domain result structur
alldom_results = cell(1,maxndom);

% plot optimizer process
plotFlag = 0;

% Monte Carlo Simulation: number of itereations
numits = 20;

% start timer for optimization
tic

% initiate plot
plot_Ddata(Ddata,1);

for ndom = minndom:maxndom % number of possible domains to be evaluated

    disp([num2str(ndom),' domains'])
    % initial starting values for the optimizer

      % let Ea be the same value for all domain searches
      Ea = 65;

      % Let LnD0aa range between 17 and 3
      if ndom == 1; lnD0aa = 15;
      elseif ndom > 1; lnD0aa = fliplr((3:(17-3)/(ndom-1):17)); end % set random lnD0aa values between 17 and 3
      
      if ndom == 1; Fx = 1;
      else 
            randomFX = sort(rand(1,ndom),'descend');
            Fx = randomFX/sum(randomFX);
      end
    
    % structure starting values
    X0 = [Mtot (Ea) lnD0aa (Fx)];
    minx = [min(Mtotlim) min(Ealim) ones(1,ndom).*min(lnD0aalim) ones(1,ndom).*min(Fxlim)];
    maxx = [max(Mtotlim) max(Ealim) ones(1,ndom).*max(lnD0aalim) ones(1,ndom).*max(Fxlim)];
        
    % get the initial X0 misfit
    M0 = objective_MDD_model(X0,Ddata,0,init_Diff);
    
    % optimizer linear constraints
    % Constraint such that D0 is kept in order and sum of Fx = 1
    if ndom > 1    
        Acon = [zeros(ndom-1,1) zeros(ndom-1,1) zeros(ndom-1,ndom)  zeros(ndom-1,ndom)];     
        Aeqcon = [0 0 zeros(1,ndom) ones(1,ndom)];
        for a = 1:ndom-1
            Acon(a,2+a) = -1;
            Acon(a,2+(a+1)) = 1;
        end
        bcon = [zeros(ndom-1,1)];
        beqcon = 1;
    else
        Acon = [];
        bcon = [];
        Aeqcon = [];
        beqcon = [];
        
    end
    
    % set up empty cells
    results_optx = zeros(numits,length(X0));
    results_fval = zeros(numits,1);
    results_exitflag = zeros(numits,1);

    for a = 1:numits
           
    % use 'for' when plotflag is turned on, or for smaller numits
    % use 'parfor' when plotflag is turned off and for large number of iterations
         
        % Display iterations
        if numits <= 100 && numits > 1
            disp([ num2str(a),' iterations out of ',num2str(numits)])
        elseif numits <= 1000
            for b = 10:10:numits
                if a == b
                    disp([ num2str(a),' iterations out of ',num2str(numits)])
                end
            end
        elseif numits <= 10000
            for b = 100:100:numits
                if a == b
                    disp([ num2str(a),' iterations out of ',num2str(numits)])
                end
            end
        end
       
        % % Generate random normal distributed diffusion concentration
        %    
        %       randomDdata.M = normrnd(randomDdata.M,randomDdata.delM);

        %       % At no point can you have a negative release of gas, therefore
        %       % Find all negative value and set them to zero
        %       randomDdata.M(find(randomDdata.M <= 0)) = 1e-20;
        
        % Do NOT Generate random normal distributed diffusion concentration
        randomDdata = Ddata;

        if plotFlag == 1
            % clear previous plots
            delete(findobj('tag','data'))
            delete(findobj('tag','Mi'))
            set(findobj('tag','misfit'),'color',[1 0.85 0.85]);
            % Initialize plot
            plot_Ddata(randomDdata,plotFlag);
        end
    

        % GlobalSearch
        % Documentation exist here
        % https://www.mathworks.com/help/gads/how-globalsearch-and-multistart-work.html#bsc9eec
        
        options = optimoptions('fmincon');
        options = optimoptions(@fmincon,'display','off','Algorithm','sqp','StepTolerance',1e-12,'Tolfun',1e-12,'MaxIterations',1000*length(X0),'MaxFunctionEvaluations',1000*length(X0));

        problem = createOptimProblem('fmincon','x0',X0,...
            'objective',@(X) objective_MDD_model(X,randomDdata,plotFlag,init_Diff), ...
            'Aineq',Acon,'bineq',bcon,'Aeq',Aeqcon,'beq',beqcon,'lb',minx,'ub',maxx,'options',options,'nonlcon',@(X) nonlcon_MDD(X,randomDdata));
                   
        [optx,fval,exitflag,outpt,allmins] = run(ms,problem,nms*length(X0)); % MultiStart

    
        % Other optimization options
        %[optx,fval,exitflag] = run(gs,problem);
        %[optx,fval,exitflag,output,grad] = fmincon(@(X) objective_MDD_model(X,randomDdata,plotFlag),X0,[],[],[],[],minx,maxx,[],opts);
        %[optx,fval,exitflag,output,grad] = fmincon(@(X) objective_MDD_model(X,randomDdata,plotFlag),X0,Acon,bcon,Aeqcon,beqcon,minx,maxx,@(X) nonlcon_MDD(X,randomDdata),opts);
    
        % Store optimizer results for each MC iteration
        results_optx(a,:) = optx;
        results_fval(a,:) = fval;
        results_exitflag(a,:) = exitflag;
    
        if plotFlag == 1
            pause(3) % this allows for you to examine the results before another MC iteration starts
        end
    
        
    end % end optimizer loop
    
    alldom_results{ndom}.optx = results_optx;
    alldom_results{ndom}.fval = results_fval;
    alldom_results{ndom}.exitflag = results_exitflag;

%%% spit out results
        
            good = find(alldom_results{ndom}.exitflag == 1);
            okay = find(alldom_results{ndom}.exitflag == 2);
            bad = find(not(ismember(alldom_results{ndom}.exitflag,[1 2])));
            
            pass = [1 2]; % ExitFlag that are allowed to pass
            
            ok = find(ismember(alldom_results{ndom}.exitflag,pass));
            
            % unpack results
           
            alldom_results{ndom}.misfit = alldom_results{ndom}.fval(ok);
            alldom_results{ndom}.Mtot = alldom_results{ndom}.optx(ok,1);
            alldom_results{ndom}.Ea = alldom_results{ndom}.optx(ok,2);
            alldom_results{ndom}.lnD0aa = alldom_results{ndom}.optx(ok,3:(ndom+2));
            alldom_results{ndom}.Fx = alldom_results{ndom}.optx(ok,ndom+3:end);
            alldom_results{ndom}.imin = find(alldom_results{ndom}.fval == min(alldom_results{ndom}.fval));
        
            % Extract best fit from each domain
        
            i = alldom_results{ndom}.imin;
            if length(i) > 1; i = i(1); end % only if more than one min is found
        
            dom_results{ndom}.fval = alldom_results{ndom}.fval(i);
            dom_results{ndom}.optx = alldom_results{ndom}.optx(i,:);
            dom_results{ndom}.Mtot = alldom_results{ndom}.optx(i,1);
            dom_results{ndom}.Ea = alldom_results{ndom}.optx(i,2);
            dom_results{ndom}.lnD0aa = alldom_results{ndom}.optx(i,3:(ndom+2));
            dom_results{ndom}.Fx = alldom_results{ndom}.optx(i,ndom+3:end);
        
            fvals(ndom,1) = dom_results{ndom}.fval;
            Eas(ndom,1) = dom_results{ndom}.Ea;
        
        
            % display best fit results
            
            disp(['',newline,'Optimization Results:'])
            disp(' ')
            disp([sprintf(['Misfit:' ...
                '\t%0.2f'],dom_results{ndom}.fval)])
            disp([sprintf('Activation Energy:\t%0.2f',dom_results{ndom}.Ea)])
            
            for a = 1:ndom
                disp([sprintf('lnD0aa %0.0f:\t\t%0.2f',a,dom_results{ndom}.lnD0aa(a))])
            end
           
            for a = 1:ndom
                disp([sprintf('Fraction %0.0f:\t\t%0.3f',a,dom_results{ndom}.Fx(a))])
            end
 
            % now plot best fit for ndom
            objective_MDD_model(dom_results{ndom}.optx,Ddata,1,init_Diff);


end % end domain loop

% Save MDD results

% report optimization time
time = toc/60/60

% Make a sound when optimization is done

S(1) = load('chirp');
S(2) = load('splat');
sound(S(1).y,S(1).Fs)
pause(1.6)
sound(S(2).y,S(2).Fs)

%% unpack all domain results
close all;

% initiate plot
plot_Ddata(Ddata,1);

for ndom = minndom:maxndom
    
    good = find(alldom_results{ndom}.exitflag == 1);
    okay = find(alldom_results{ndom}.exitflag == 2);
    bad = find(not(ismember(alldom_results{ndom}.exitflag,[1 2])));
    
    pass = [1 2]; % ExitFlag that are allowed to pass
    
    ok = find(ismember(alldom_results{ndom}.exitflag,pass));
    
    % unpack results
   
    alldom_results{ndom}.misfit = alldom_results{ndom}.fval(ok);
    alldom_results{ndom}.Mtot = alldom_results{ndom}.optx(ok,1);
    alldom_results{ndom}.Ea = alldom_results{ndom}.optx(ok,2);
    alldom_results{ndom}.lnD0aa = alldom_results{ndom}.optx(ok,3:(ndom+2));
    alldom_results{ndom}.Fx = alldom_results{ndom}.optx(ok,ndom+3:end);
    alldom_results{ndom}.imin = find(alldom_results{ndom}.fval == min(alldom_results{ndom}.fval));

    % Extract best fit from each domain

    i = alldom_results{ndom}.imin;

    if length(i) > 1; i = i(1); end % only if more than one min is fou
    
    dom_results{ndom}.fval = alldom_results{ndom}.fval(i);
    dom_results{ndom}.optx = alldom_results{ndom}.optx(i,:);
    dom_results{ndom}.Mtot = alldom_results{ndom}.optx(i,1);
    dom_results{ndom}.Ea = alldom_results{ndom}.optx(i,2);
    dom_results{ndom}.lnD0aa = alldom_results{ndom}.optx(i,3:(ndom+2));
    dom_results{ndom}.Fx = alldom_results{ndom}.optx(i,ndom+3:end);

    fvals(ndom,1) = dom_results{ndom}.fval;
    Eas(ndom,1) = dom_results{ndom}.Ea;


    % display best fit results
    
    disp(['',newline,'Optimization Results:'])
    disp(' ')
    disp([sprintf(['Misfit:' ...
        '\t%0.2f'],dom_results{ndom}.fval)])
    disp([sprintf('Activation Energy:\t%0.2f',dom_results{ndom}.Ea)])
    
    for a = 1:ndom
        disp([sprintf('lnD0aa %0.0f:\t\t%0.2f',a,dom_results{ndom}.lnD0aa(a))])
    end
   
    for a = 1:ndom
        disp([sprintf('Fraction %0.0f:\t\t%0.3f',a,dom_results{ndom}.Fx(a))])
    end

    % now plot best fit for ndom
    objective_MDD_model(dom_results{ndom}.optx,Ddata,1,init_Diff);
    pause(2) % use 'pause' to allow time to observe the results and continue by click any key

end

% Now evaluate optimal numbers of domain to use
format long

% % Determine if Domain misfit gets statistically better
% if numits == 1
%     dom = minndom;
%     for a = minndom+1:(maxndom)
%         if alldom_results{a}.fval < alldom_results{a-1}.fval
%             dom = a;
%         else
%             break
%         end
%     end
% 
% elseif numits >1
%     dom = minndom;
%     for a = minndom+1:(maxndom)
% 
%         if ttest2(alldom_results{a}.fval,alldom_results{(a-1)}.fval) &&  dom_results{a}.fval < dom_results{a-1}.fval
%             dom = a;  
%         else
%             break % stop if no significant difference in misfit
%         end
%     end
% end

figure(3); hold on;

for ndom = minndom:maxndom 
    plot(ones(numits,1)*ndom,alldom_results{ndom}.fval,'ob','HandleVisibility','off');
    plot(ndom,dom_results{ndom}.fval,'ob','Markerfacecolor','b','HandleVisibility','off'); 

    errorbar(ndom,mean(alldom_results{ndom}.fval),std(alldom_results{ndom}.fval),'sk','MarkerSize',8,'HandleVisibility','off'); hold on;

end

%plot(dom,dom_results{dom}.fval,'ko','markerfacecolor','g','displayname');

% legend plot
plot(NaN,NaN,'ob','displayname','All misfit values');
plot(NaN,NaN,'ob','Markerfacecolor','b','displayname','Bestfit'); 
plot(NaN,NaN,'sk','markersize',8,'displayname','mean misfit'); 

legend('boxoff');


set(gca,'YScale','log');
xlim([0 ndom])
%curtick = get(gca, 'YTick');
%set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xlabel('Number of Domains')
ylabel('Misfit value')
grid on

pause(2)

% find number of optimal domains

doms = 1:maxndom;
for a = minndom:maxndom
    fval(a) = dom_results{a}.fval; 
    df(a) = sum(Ddata.use)-length(dom_results{a}.optx);

end

% calculate reduced Chi^2 results
rfval = fval./df;


dom = find(rfval == min(rfval));

% now plot reduced misfit statistics
figure(4); hold on;

plot(doms,rfval,'.k','MarkerSize',20,'HandleVisibility','off');
plot(dom,rfval(dom),'.r','Markersize',20,'displayname','Best domain fit');

xlabel('number of domains')
ylabel('reduced misfit')

legend('boxoff');

% Display number of domain and best fit
disp([sprintf('Suggested number of domains: %0.0f',dom)])

pause(1)

objective_MDD_model(dom_results{dom}.optx,Ddata,1,init_Diff);

%% MDDdata

MDDdata = dom_results{dom};
MDDdata.D0aa = exp(MDDdata.lnD0aa);

f = fieldnames(Ddata);
 for a = 1:length(f)
    MDDdata.(f{a}) = Ddata.(f{a});
 end

MDDdata.init_Diff = init_Diff;

MDDdata = plot_MDD_forward_model(MDDdata,Ddata,1,init_Diff);

%% Ask to save MDD Results

%cd '/Users/marie/Desktop/MatLab/Diffusion'

if isfile(fullfile(cd, "./Ddata/MDDdom_results/"+Ddata.samplename+ ".mat"))

    q = questdlg('Do you want to override MDD Results?','Override MDD Results','Yes','No','No');
    switch q
        case 'Yes'
        
            % save workspace as a .mat file
            save("./Ddata/MDDdom_results/"+Ddata.samplename+ ".mat");
            % Save MDD results as a structure for plotting
            save("./Ddata/MDDdata/"+Ddata.samplename+ ".mat","MDDdata")

            disp('MDD Results sucessfully saved')
    
        case 'No'  
            disp('MDD Results not saved')
    end
else

        q = questdlg('Do you want to save MDD Results?','Save MDD Results','Yes','No','No');
    switch q
        case 'Yes'
        
            % save workspace as a .mat file
            save("./Ddata/MDDdom_results/"+Ddata.samplename+ ".mat");
            % Save MDD results as a structure for plotting
            save("./Ddata/MDDdata/"+Ddata.samplename+ ".mat","MDDdata")

            disp('MDD Results sucessfully saved')
    
        case 'No'  
            disp('MDD Results not saved')
    end
end


