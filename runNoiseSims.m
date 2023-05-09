function [S] = runNoiseSims(P)
% Runs multiple simulations with the same input profile, with option for added white
% noise to the input.
%
% Inputs:
% P - parameter structure with fields indicating simulation parameters,
% the neuron model type to use, and input current parameters. Fields listed below
% with example values.
%
% Outputs:
% S - Simulation data structure.

% % Example input params
% P.t0 = 8; % start time of input. Make sure this is long enough to allow V and U 
% % to reach steady states before input arrives.
% P.A = 500; % Slope of the tent input.  
% P.beta = 1000; % Height of the tent input peak.
% P.sigma = 10; % Noise amplitude
% 
% % Simulation params struct
% P.model = 'Meng VU';
% P.input = 'tent';
% P.numSims = 1; % Number of simulations to run.
% P.numT = 2; % Number of input periods to run per simulation.
% P.dt = 0.001;
% P.simLength = P.t0 + P.numT*4*P.beta/P.A;
% P.tarray = [0 : P.dt : P.simLength]';
% P.numSteps = length(P.tarray);
% P.spike_thr = 0; % Spike detection voltage threshold;
% P.reset_thr = -40; % Threshold to reset spike detection;
% P.V0 = -60; % voltage initial condition in mV
% P.U0 = 0; % adaptation variable initial condition


switch P.model

    case {'VU phasic', 'VU tonic'}
        % From Meng, Huguet, Rinzel 2012
        w0 = 0.511; h0 = 0.445; r0 = 0.147; z0 = 0.662; n0 = 0.0077; p0 = 0.0011;
        a = 0.9; b = (a - w0)/h0;
        C=12; gna=1000; gkht=150; gklt=200; gh=20; glk=2;
        ENa=55; EK=-70; Eh=-43; Elk=-65;
        hhalf=-65; whalf=-48; sigh=6; sigw=6;
        
        hinf = @(V) 1.0/(1.0+exp((V-hhalf)/sigh));
        winf = @(V) (1.0/(1.0+exp(-(V-whalf)/sigw)))^(0.25);
        uinf = @(V) b*(hinf(V)+b*(a-winf(V)))/(a*(1+b^2));
        minf = @(V) 1.0/(1.0+exp(-(V+38)/7));
        
        tauw = @(V) 100/(6*exp((V+60)/6)+16*exp(-(V+60)/45))+1.5;
        tauh = @(V) 100/(7*exp((V+60)/11)+10*exp(-(V+60)/25))+0.6;
        tau = @(V) min(tauw(V),tauh(V));

    otherwise
        error("Invalid neuron model type provided.")
end


% Arrays to hold data across simulations
V_mat = zeros(P.numSteps, P.numSims);
U_mat = zeros(P.numSteps, P.numSims);
input_mat = zeros(P.numSteps, P.numSims);
spikes_cell = cell(1,P.numSims); % Times at which spikes occur

input_noNoise = zeros(P.numSteps,1);
T_t0 = P.t0; % For the first period that is extra long to allow ICs to settle.

% Run a loop to get the input profile without noise
for i = 1:(P.numSteps - 1)
    % To run multiple periods of input I need to update the second arg in
    % input() to the current time once a new period starts.
    
    T_newt0 = P.T_startTimes(find(P.dt*i > P.T_startTimes, 1, 'last')); % Finds last element
    % in T_startTimes that is less than current time dt*i.

    if isempty(T_newt0) == 0 && T_newt0 ~= T_t0
        T_t0 = T_newt0; % Updates start time of period
    end

    % adds phase shifted second input profile.
    if P.isSecondInput
        input_noNoise(i) = P.input(P.dt*i,T_t0) + P.input(P.dt*i, T_t0+P.inputShift);
    elseif ~P.isSecondInput % single input profile
        input_noNoise(i) = P.input(P.dt*i,T_t0);
    else
        error("Invalid P.isSecondInput. Only 0 and 1 are accepted.  ")
    end

    

end


% Run actual simulations
for sim = 1:P.numSims

    % Arrays to hold simulation history
    V_sim = zeros(P.numSteps,1);
    U_sim = zeros(P.numSteps,1);
    input_sim = zeros(P.numSteps,1);
    spikes_sim = [];
    
    V_sim(1) = P.V0; U_sim(1) = P.U0; % Set ICs
    spike_boo = 0; % Boolean for if spike has occured.

    % Run simulation
    for i = 1:(P.numSteps - 1)

        V = V_sim(i); U = U_sim(i);

        I = input_noNoise(i) + (i*P.dt > P.t0)*P.sigma*randn()/sqrt(P.dt); % Input current with noise
        input_sim(i) = I;
    
        % Detects a spike if voltage above threshold, spike has not just
        % occured, and input has started
        if V >= P.spike_thr && spike_boo == 0 && P.dt*i >= P.t0
            spikes_sim(end+1) = P.dt*i; % Records time (in ms) of spike
            spike_boo = 1;
        elseif V <= P.reset_thr && spike_boo == 1
            spike_boo = 0; % Resets so ready to detect another spike.
        end
    
        switch P.model
            case 'VU phasic'
                dV = P.dt*((2*(-gna*minf(V)^3*(a*U/b)*(V-ENa)-gklt*a^4*(1-U)^4*z0*(V-EK)-gkht*(0.85*n0^2+0.15*p0)*...
                    (V-EK)-glk*(V-Elk)-gh*r0*(V-Eh))+I)/C);
                dU = P.dt*(3*(uinf(V)-U)/tau(V));
            case 'VU tonic' % U is frozen to IC value in the gklt term.
                dV = P.dt*((2*(-gna*minf(V)^3*(a*U/b)*(V-ENa)-gklt*a^4*(1-P.U0)^4*z0*(V-EK)-gkht*(0.85*n0^2+0.15*p0)*...
                    (V-EK)-glk*(V-Elk)-gh*r0*(V-Eh))+I)/C);
                dU = P.dt*(3*(uinf(V)-U)/tau(V));
        end
    
        V_sim(i+1) = V + dV; U_sim(i+1) = U + dU;
    
    
        
    end
    
    V_mat(:,sim) = V_sim; U_mat(:,sim) = U_sim; % Stores simulation.
    input_mat(:,sim) = input_sim;
    spikes_cell(sim) = {spikes_sim};

end

S.V_mat = V_mat; S.U_mat = U_mat;
S.input_mat = input_mat;
S.spikes_cell = spikes_cell;
S.input_noNoise = input_noNoise;

nonempty_inds = ~cellfun('isempty',spikes_cell); % Identify indices of sims 
% that have a spike.
first_spikes = cellfun(@(c) c(1), spikes_cell(nonempty_inds)); % Extract 
% only the first spike time from all cells with spikes.
S.first_spikes = first_spikes;


end