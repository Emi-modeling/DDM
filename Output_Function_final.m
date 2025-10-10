%rt = [0.45, 0.60, 0.38, 0.52, 0.70];       % vector of RTs in seconds
%acc = [1, 0, 1, 1, 0] ;      % vector of accuracy (0 or 1)
%Instr = {'F', 'F', 'I+F', 'F', 'I+F'};    % cell array or categorical of instruction conditions
%initial values for parameters = priors
%pars  : parameter vector [v1 eta1 a1 sz1 Ter1 st01 z1 v2 eta2 a2 sz2 Ter2 st02 z2 ]
%pars = [0.3 0.1 1.2 0.1 0.2 0.05 0.5,  0.5 0.1 1.0 0.1 0.25 0.05 0.5];


% these example data are random and the trial number is too small to make
% any interpretation of parameters or fit!

% alternatively simulate data based on DDM parameters

nTrialsPerCond = 200;
conditions = {'F', 'I+F'};
nCond = numel(conditions);

% Parameter per condition (v, eta, a, sz, Ter, st0, z)
params = [
    0.3, 0, 1.2, 0.1, 0.3, 0.05, 0.5;   % F
    0.5, 0, 1.2, 0.1, 0.3, 0.05, 0.5    % I+F
];

% Initializing vectors 
rt_all = [];
acc_all = [];
instr_all = {};

% function for simulation
for c = 1:nCond
    par = params(c, :);
    v = par(1); a = par(3); sz = par(4)*a; Ter = par(5); st0 = par(6); zrel = par(7);
    z = zrel * a;

    for t = 1:nTrialsPerCond
        % Drift + Noise
        drift = v;
        t_eff = 0;
        x = z;
        dt = 0.001;
        noise = 1;
        maxTime = 5;  % timeout in seconds
        upper = a;
        lower = 0;
        
        while true
            dx = drift * dt + noise * sqrt(dt) * randn;
            x = x + dx;
            t_eff = t_eff + dt;
            if x >= upper
                rt = Ter + st0 * rand + t_eff;
                acc = 1;
                break
            elseif x <= lower
                rt = Ter + st0 * rand + t_eff;
                acc = 0;
                break
            elseif t_eff > maxTime
                rt = Ter + st0 * rand + t_eff;
                acc = rand > 0.5;
                break
            end
        end
        
        rt_all(end+1,1) = rt;
        acc_all(end+1,1) = acc;
        instr_all{end+1,1} = conditions{c};
    end
end

% show results
fprintf('Gesamttrials: %d\n', numel(rt_all));
fprintf('Mean RT (F): %.3f\n', mean(rt_all(strcmp(instr_all, 'F'))));
fprintf('Accuracy (F): %.2f\n', mean(acc_all(strcmp(instr_all, 'F'))));
fprintf('Mean RT (I+F): %.3f\n', mean(rt_all(strcmp(instr_all, 'I+F'))));
fprintf('Accuracy (I+F): %.2f\n', mean(acc_all(strcmp(instr_all, 'I+F'))));

% save results
save('simDDMdata.mat', 'rt_all', 'acc_all', 'instr_all');



% Objective function for PSO which has to be minimized
objFun = @(params) negLL_DDM_final(params, rt_all, acc_all, instr_all); 




