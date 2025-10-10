function negLL = negLL_DDM_multi(pars, rt, acc, Instr)
    % negLL_DDM_multi - Full DDM Negative Log-Likelihood with all parameters
    % Inputs:
    %   pars : [v, eta, a, sz, Ter, st0, z]_1, ..., _nCond
    %   rt   : reaction times
    %   acc  : accuracy (0/1)
    %   Instr: condition labels
    %
    % Output:
    %   negLL: total negative log-likelihood
    
    response = acc + 1;  % DDM expects 1 = lower (error), 2 = upper (correct)
    uniqueInstr = unique(Instr);
    nCond = numel(uniqueInstr);
    nParsPerCond = length(pars) / nCond;
    negLL = 0;

    for c = 1:nCond
        idx = (c-1)*nParsPerCond + (1:nParsPerCond);
        condPars = pars(idx);
        condMask = strcmp(Instr, uniqueInstr{c});
        negLL = negLL + negLL_DDM_single_full(condPars, rt(condMask), response(condMask));
    end
end


function negLL = negLL_DDM_single_full(pars, rt, response)
    % Unpack parameters
    v    = pars(1);
    eta  = pars(2);
    a    = pars(3);
    sz   = pars(4) * a;
    Ter  = pars(5);
    st0  = pars(6);
    zrel = pars(7);
    z    = zrel * a;

    % Remove non-decision time
    nTrials = length(rt);
    rt_adj = rt - Ter;
    rt_adj(rt_adj <= 0) = eps;

    densities = zeros(size(rt_adj));
    
    % Sampling steps for marginalization
    nSamples = 5; % tradeoff between accuracy and speed

    for i = 1:nTrials
        total = 0;
        for s = 1:nSamples
            % Sample z from uniform distribution (start point variability)
            if sz > 0
                z_sample = z + (rand - 0.5) * sz;
                z_sample = min(max(z_sample, 0), a);
            else
                z_sample = z;
            end

            % Sample Ter from uniform distribution (non-decision variability)
            if st0 > 0
                Ter_sample = Ter + rand * st0;
            else
                Ter_sample = Ter;
            end

            t_sample = rt(i) - Ter_sample;
            if t_sample <= 0
                continue;
            end

            % Sample drift rate from normal distribution (eta)
            if eta > 0
                v_sample = v + eta * randn;
            else
                v_sample = v;
            end

            % DDM PDF (approximate)
            p = ddiffusion_pdf(t_sample, response(i), a, v_sample, z_sample);
            total = total + p;
        end

        densities(i) = total / nSamples;
    end

    % Contaminant (uniform) noise
    rt_min = min(rt_adj);
    rt_max = max(rt_adj);
    contaminantDensity = 1 / (rt_max - rt_min);
    finalDensity = 0.05 * contaminantDensity + 0.95 * densities;

    % Prevent log(0)
    finalDensity(finalDensity < 1e-10) = 1e-10;

    % Negative log-likelihood
    negLL = -sum(log(finalDensity));
end


function p = ddiffusion_pdf(t, response, a, v, z)
    % Approximation of first-passage time density for DDM
    s = 1; % diffusion constant
    nTerms = 30;

    % Reflect if response is lower boundary
    if response == 1
        v = -v;
        z = a - z;
    end

    if t <= 0
        p = 0;
        return;
    end

    w = z / a;
    t_scaled = (a^2 / s^2) * t;
    sumTerm = 0;

    for k = -nTerms:nTerms
        m = 2 * k + w;
        term = m * exp(-((m^2) * pi^2 * t_scaled) / 2) * sin(m * pi);
        sumTerm = sumTerm + term;
    end

    p = (pi / a^2) * sumTerm;

    if isnan(p) || p < 0
        p = 0;
    end
end
