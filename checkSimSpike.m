function [] = checkSimSpike(P)
% Runs the baseline simulation and raises an error if a spike occurs.
%
% Inputs:
% P - parameter struct with simulation specifications. More details about P
% provided in runNoiseSims.m.
%
% Outputs: None
% 


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



% Arrays to hold baseline simulation history
V_sim = zeros(P.numSteps,1);
U_sim = zeros(P.numSteps,1);
input_sim = zeros(P.numSteps,1);

V_sim(1) = P.V0; U_sim(1) = P.U0; % Set ICs
% Run baseline simulation to make sure there is no spiking due to baseline
% current.
for i = 1:(P.numSteps - 1)
    V = V_sim(i); U = U_sim(i);

    I = P.input(P.dt*i,P.t0); % Input current
    input_sim(i) = I;

    % Throws error if spike occurs. Spikes should not occur during baseline
    % simulation because then no kick is needed and the dynamical threshold
    % cannot be measured.
    if V >= P.spike_thr
        error("Spike occurred during baseline simulation using A = %d, beta = %d.", P.A, P.beta)
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

end