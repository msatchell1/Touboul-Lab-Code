function [DT] = getDynThr(P)
% Measures the dynamical threshold for a given set of parameters P. Note
% this function does not currently support multiple periods of input.
%
% Inputs:
% P - parameter structure with fields indicating simulation parameters,
% the neuron model type to use, and input current parameters. Fields listed below
% with example values.
%
% Outputs:
% DT - Dynamical threshold structure.
% 
% % Example input current params
% P.t0 = 8; % start time of input. Make sure this is long enough to allow V and U 
% % to reach steady states before input arrives.
% P.A = 500; % Slope of the tent input.  
% P.beta = 700; % Height of the tent input peak.
% 
% % Example simulation params
% P.model = 'Meng VU';
% P.input = 'tent';
% P.dt = 0.001;
% P.simLength = P.t0 + 4*P.beta/P.A;
% P.tarray = [0 : P.dt : P.simLength]';
% P.numSteps = length(P.tarray);
% P.spike_thr = 0; % Spike detection voltage threshold;
% P.reset_thr = -40; % Threshold to reset spike detection;
% P.V0 = -60; % voltage initial condition in mV
% P.U0 = 0; % adaptation variable initial condition
% 
% % Example kicking params
% P.kickDensity = 15; % How many kicks to have per milisecond.
% % Evenly spaced kick times of P.kickDensity density. Note kickTimes can't
% % start at 0 because I get an indexing by 0 error later.
% P.kickTimes = linspace(P.dt, P.simLength, P.kickDensity*(P.simLength-P.dt));
% P.kickIncrement = 0.1; % Size to increase kick amplitude each iteration.



switch P.model

    case {'VU phasic','VU tonic'}
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


% Determines whether to use 1 or 2 input profiles.
if P.isSecondInput % adds phase shifted second input profile.
    DT_input = @(i,T_t0) P.input(P.dt*i,T_t0) + P.input(P.dt*i, T_t0+P.inputShift);
elseif ~P.isSecondInput % single input profile
    DT_input = @(i,T_t0) P.input(P.dt*i,T_t0);
else
    error("Invalid P.isSecondInput. Only 0 and 1 are accepted.  ")
end



% Arrays to hold baseline simulation history
V_sim = zeros(P.numSteps,1);
U_sim = zeros(P.numSteps,1);
input_sim = zeros(P.numSteps,1);

V_sim(1) = P.V0; U_sim(1) = P.U0; % Set ICs

% Run baseline simulation
for i = 1:(P.numSteps - 1)
    V = V_sim(i); U = U_sim(i);

    I = DT_input(i,P.t0); % gets input current
    input_sim(i) = I;

    % Throws error if spike occurs. Spikes should not occur during baseline
    % simulation because then no kick is needed and the dynamical threshold
    % cannot be measured.
    if V >= P.spike_thr
        error("Spike occurred during baseline simulation.")
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




% Arrays to hold data across simulations
spike_times = zeros(1, length(P.kickIndices)); % Times at which spikes occur
kick_sizes = zeros(1, length(P.kickIndices)); % Kick sizes that elicit a spike in mV.


% Actual dynamical threshold sims
for j = 1:length(P.kickIndices)
    
    if j == 1
        kick = 0; % Size of kick in mV
    else % to try and make the code more efficient, I am trying not to always
        % start the kicksize search from zero. It seems that most DTs don't
        % have slopes over 40, so I should be safe assuming a floor of
        % 50/KD. I should only ever run into a problem if I calculate an
        % extremely sharp input profile with a very low KD.
        kick = max( kick_sizes(j-1) - 1 - 50/P.kickDensity, 0); % max(,0) ensures no 
        % negative kick voltages.
    end

    spike_boo = 0;
    kick_i = P.kickIndices(j); % Index to kick at
    V_IC = V_sim(kick_i); U_IC = U_sim(kick_i); % Initial conditions

    while ~spike_boo

        kick = kick + P.kickIncrement; % Increase kick size
        V = V_IC + kick; % Apply kick to voltage
        U = U_IC;
        i = kick_i; % Starting index for simulation
        
        while i <= (P.numSteps)

            if V >= P.spike_thr % If spike occurs

                % To get the most accurate spike time I will interpolate
                % linearly based on the two voltages on either side of the
                % spike threshold.
                interp_i = interp1([V-dV,V],[i-1,i],P.spike_thr,"linear");

                kick_sizes(1,j) = kick;
                spike_times(1,j) = interp_i*P.dt; % Time of spike

                spike_boo = 1; % Stops kick loop
                i = P.numSteps; % Stops simulation 
            end

            I = DT_input(i,P.t0);
            
%             % I could put switch-case here to allow for more models... But
%             % I am afraid this would affect run time.
%             dV = P.dt*((2*(-gna*minf(V)^3*(a*U/b)*(V-ENa)-gklt*a^4*(1-U)^4*z0*(V-EK)-gkht*(0.85*n0^2+0.15*p0)*...
%                 (V-EK)-glk*(V-Elk)-gh*r0*(V-Eh))+I)/C);
%             dU = P.dt*(3*(uinf(V)-U)/tau(V));

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

            V = V + dV; U = U + dU;
            i = i+1;
        end
       

    end

    fprintf("Kick time %d out of %d complete \n", j, length(P.kickTimes))

end

% Dynamical threshold data structure
DT.kickDensity = P.kickDensity;
DT.kickTimes = P.kickTimes;
DT.kickIncrement = P.kickIncrement;
DT.kickSizes = kick_sizes;
DT.spikeTimes = spike_times;
DT.input = input_sim;
DT.spike_lags = spike_times - P.kickTimes;



end