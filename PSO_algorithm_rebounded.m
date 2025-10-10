% Main program for a Particle Swarm Algorithm like (Genetic algorithm like)
% Variation realized by Paolo Mercorelli, Constanze Fuchs and Emelie Lenze
% from the original algorithm realized by Paolo Mercorelli
% ABB-Research Center Heidelberg in the context of non convex optimization algorithms for 
%Pantograph Control: Patent number DE20021058921 20021217  

% initialization algorithm
tic % tic and toc measures time between two commands 
%clc % clears command window 
%clear all
%close all

rng default % randomly moving particles/ direction of vectors, uniform distribution

format long



% initializing instructions and conditions
uniqueInstr = unique(instr_all);
nCond = numel(uniqueInstr);

% initializing boundaries for the parameters 
% [v, eta, a, sz, Ter, st0, z]
alpha = 0.1;
LB_single = [-5,  0,   0.3,  0,    0.1,  0,   0];
UB_single = [ 5,  0.5, 3.0,  0.5,  0.6,  0.3, 1];
% initializing weights vector to rebound boundaries 10% of their range 
Weights = [abs(LB_single-UB_single)*alpha, abs(LB_single-UB_single)*alpha];
%%
LB = repmat(LB_single, 1, nCond);
UB = repmat(UB_single, 1, nCond);

% transformed boundaries
%for i = 1:length(LB)
    %LB(i) = param_transform(LB(i), transform_all{i});
    %UB(i) = param_transform(UB(i), transform_all{i});
%end


% Parameters of Particle Swarm Algorithm

m= 14;  %number of variables to be optimized, 2 conditions with 7 parameters each

n=100; % population size, using more particles makes solution better

wmax=0.9; % inertia weight

wmin=0.4; % inertia weight

k1=2;  %acceleration factor, can also be variable 

k2=2; %acceleration factor

intermeadian_iteration=1000; % how to stop

% Main Particle Swarm Optimization Algorithm

maxiteration=1000; % maximal number of iterations

maxrunloops=1;     %maximal number of loops



global xdata ydata

%load 

for run=1:maxrunloops
run
for i=1:n
    
    for j=1:m
        
        x0(i,j)=(LB(j)+rand()*(UB(j)-LB(j))); %fill up between lower and upper bound the position of particles= initial start position of swarm
        
        %x0(i,j)=round((LB(j)+UB(j))/2+((UB(j)-LB(j)))*rand());
    end
end
end

x=x0;  % initial matrix of population position

v=0.1*x0; % initial matrix population velocity

for i=1:n
    
   f0(i,1)=objFun(x0(i,:));  % calling output function, all values are in a column, using randomly definded x0 
end
    
[fmin0,index0]=min(f0);   % fmin= find the minimum of this column, index0= where it is, store it in element min(f0)

%pause
pbest=x0;    % initial pbest (position best), matrix

gbest=x0(index0,:);  % initial gbest (global best = best solution for initial postition of particles) 

% PSO Main program genetic algorithm

iteration=1;

tollerance=1;

while iteration<=maxiteration;  % && tollerance>10^-6

    w=wmax-(wmax-wmin)*iteration/maxiteration;  %update weights = during iteration, we try to converge at one point =explotation 
    
    % velocity parameters update
    
    %k1=wmax-(wmax-wmin)*iteration/maxiteration;
    
   % k2=wmax-(wmax-wmin)*iteration/maxiteration;
    
    for i=1:n
%     
    for j=1:m
%         

         v(i,j)=w*v(i,j)+k1*rand()*(pbest(i,j)-x(i,j))+k2*rand()*(gbest(1,j)-x(i,j));
         % matrix reduces every iteration, 
         % vectors change every iteration, 
         % adjust pbest = all posotion of all particles 
         % adjsut gbest = bird that tries to orientate the swarm because of best
         % intial position 
     end
    end 

   for i=1:n
    
    for j=1:m
        
        x(i,j)=x(i,j)+v(i,j);
    end
   end 
   
   % check boundaries
   
   for i=1:n;
      
       for j=1:m;
           
           if x(i,j)<LB(j)
               x(i,j)=LB(j);
               v(i,j)=Weights(1,j)*rand();
               x(i,j)=x(i,j)+v(i,j);
           elseif x(i,j)>UB(j)
            x(i,j)=UB(j);  
            v(i,j)=-Weights(1,j)*rand();
            x(i,j)=x(i,j)+v(i,j);
       end
       end
   end
   
   % evaluation of the cost function
   
   for i=1:n
       
       f(i,1)=objFun(x(i,:));  % all values are ordered in a column
   end
   
   % update pbest and const function
   
   for i=1:n
           
           if f(i,1)<f0(i,1)  % if new positon is less as initial position
               pbest(i,:)=x(i,:); % pbest is replaced, for each particle with all components (x1, x2, x3)
          f0(i,1)=f(i,1);   % f is replaced
       end
   end
       
   [fmin,index]=min(f0); % find the best particle with value and position
   
   bestfmin(iteration,run)=fmin; % store the best cost fuction at the iterationa and run index
   
   bestiteration(run)=iteration; % store the iteration count in which we found the best
   
    % update gbest and best const function
    
   %for j=1:n;
           
           if fmin<fmin0
               gbest=pbest(index,:); %replace best of the best with index = position
               %fmins=fmin;
          fmin0=fmin;   
       end
   %end
    
   % check the tollerance
   
   if iteration>intermeadian_iteration
       
       tollerance=abs(bestfmin(intermeadian_iteration,run)-fmin0); % Adjusting tollernace 
    %cond1=bestfmin(intermeadian_iteration,run);
   
  
   
  
   
  if fmin0>0 && (tollerance)>0 && (tollerance)<1e-3 
% 
        gbest
%        
       
           
          

bestPars = gbest;  % Best-fit DDM parameters
bestNegLL = objFun(bestPars);

        return
%    %end   
    end
     end
     
      iteration=iteration+1;
 end 

%k = numel(bestPars);  % Number of parameters
%n = numel(rt_all);    % Number of trials
%AIC = 2*k + 2*bestNegLL;
%BIC = k*log(n) + 2*bestNegLL;

bestPars = gbest;  % Best-fit DDM parameters
bestNegLL = objFun(bestPars);  

 gbest 
 bestNegLL
 %AIC
 %BIC

toc

%% Visual Fit Check

% Parameters
nSimTrials = 200;  % Number of simulated trials per condition
simRT = [];
simAcc = [];
simCond = {};

uniqueInstr = unique(instr_all);
nCond = numel(uniqueInstr);
nParsPerCond = numel(bestPars) / nCond;

for c = 1:nCond
    % Extract best-fitting parameters for this condition
    idx = (c-1)*nParsPerCond + (1:nParsPerCond);
    par = bestPars(idx);
    
    v = par(1);
    eta = par(2);
    a = par(3);
    sz = par(4) * a;
    Ter = par(5);
    st0 = par(6);
    z = par(7) * a;

    for t = 1:nSimTrials
        % Drift sample
        if eta > 0
            drift = v + eta * randn;
        else
            drift = v;
        end

        t_eff = 0;
        x = z;
        dt = 0.001;
        noise = 1;
        maxTime = 5;
        upper = a;
        lower = 0;
        
        while true
            dx = drift * dt + noise * sqrt(dt) * randn;
            x = x + dx;
            t_eff = t_eff + dt;
            if x >= upper
                rt = Ter + st0 * rand + t_eff;
                acc = 1;
                break;
            elseif x <= lower
                rt = Ter + st0 * rand + t_eff;
                acc = 0;
                break;
            elseif t_eff > maxTime
                rt = Ter + st0 * rand + t_eff;
                acc = rand > 0.5;
                break;
            end
        end

        simRT(end+1,1) = rt;
        simAcc(end+1,1) = acc;
        simCond{end+1,1} = uniqueInstr{c};
    end
end

%% Plotting Results

figure;
for c = 1:nCond
    cond = uniqueInstr{c};
    
    % Extract empirical
    realRT = rt_all(strcmp(instr_all, cond));
    realAcc = acc_all(strcmp(instr_all, cond));
    
    % Extract simulated
    simRT_cond = simRT(strcmp(simCond, cond));
    simAcc_cond = simAcc(strcmp(simCond, cond));
    
    % RT Histogram
    subplot(nCond, 2, (c-1)*2 + 1);
    histogram(realRT, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Real');
    hold on;
    histogram(simRT_cond, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Simulated');
    title(['RT Distribution - ' cond]);
    xlabel('RT (s)');
    ylabel('Probability');
    legend;
    box on;

    % Accuracy Bar Plot
    subplot(nCond, 2, (c-1)*2 + 2);
    bar([1, 2], [mean(realAcc), mean(simAcc_cond)]);
    set(gca, 'XTickLabel', {'Real', 'Simulated'});
    ylim([0, 1]);
    ylabel('Accuracy');
    title(['Accuracy - ' cond]);
    box on;
end



