% Main script to execute simulation of a single neuron. Goals of this
% script were originally to:
% 1) Simulate a neuron model of my choice.
% 2) Provide different inputs of my choice to the simulation.
% 3) Run this simulation for long enough (or rerun the simulation many
% times) with white noise to get a histogram of spike times with respect to the period of
% input oscillation. 
% 4) Measure the dynamical threshold for a given input by adding kicks of 
% varying amplitude. For a given input (non-spike inducing, I believe), at
% each time point over one input period I will
% apply increasing kick amplitudes. After a kick is applied, I run the
% simulation out for at least one (maybe more) period to see if a spike
% occurs at any point. Once I have found the
% amplitude that induces a spike, I record that amplitude
% as the "dynamical threshold" for the kicking time point in the period, and then
% proceed on to the next time step. Repeating this process produces a
% dynamical threshold curve that shows how far the neuron is from spiking
% at each point in the input. For the V-U Meng et al. 2012 model and other
% adaptive models, the times closest to spiking are often not the times of
% maximum current amplitude, but rather those where the slope of the input
% is greatest. 
% 5) Plot all of these together to show the firing of the neuron with
% noise, the input profile, and the dynamical threshold. 

%...
% Additional sections have been added for other simulations and analyses.

%% runNoiseSims

% Input params
P.t0 = 1; % start time of input. Make sure this is long enough to allow V and U 
% to reach steady states before input arrives.
P.inputShift = 0; % The shift in time between the start of the first input and 
% the start of the second input in a single period, when there are two
% inputs being provided. In ms.
P.isSecondInput = 1; % boolean for second input current. Must only be 0 or 1.
P.A = 400; % Slope of the tent input.  
P.beta = 200; % Height of the input peak for both sin and tent.
P.sigma = 96; % Noise amplitude
P.sin_width = 10; % width of single positive sin peak.

% Simulation params struct
P.model = 'VU tonic'; % Options: 'VU phasic' or 'VU tonic'
P.input_type = 'rectsin'; % Options: 'tent' or 'rectsin'

switch P.input_type
    case 'tent'
        % Input as a tent function
        tent = @(t,t0,A,beta) (t<(t0+beta/A))*A*(t-t0) + (t>=(t0+beta/A))*(-A*(t-(t0+2*beta/A)));
        % The tent function drops down to negative instead of stopping at zero, so
        % we always need to take the maximum value between the tent function and 0
        % to prevent negative inputs.
        P.input = @(t,t0) (t>t0)*(t<=(t0+2*P.beta/P.A))*tent(t,t0,P.A,P.beta);
        P.Tlength = 2*P.beta/P.A + P.isSecondInput*(P.inputShift + 2*P.beta/P.A); % Length of a single period (ms).
    case 'rectsin'
        % input as a rectified sine wave. Width represents the width of the
        % single positive half of the wave, and so is half the length of
        % what the total period would be.
        sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
        P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
        P.Tlength = P.sin_width + P.isSecondInput*(P.inputShift + P.sin_width); % Length of a single period (ms).
    otherwise
        error("Invalid input type provided.")
end

P.numT = 1; % Number of input periods to run per simulation.
P.simLength = P.t0 + P.numT*P.Tlength;
% Period start times. Each period begins a new set of input(s).
P.T_startTimes = P.t0 : P.Tlength : P.simLength - P.Tlength;
P.numSims = 1000; % Number of simulations to run.
P.dt = 0.001;
P.tarray = [0 : P.dt : P.simLength]';
P.numSteps = length(P.tarray);
P.spike_thr = 0; % Spike detection voltage threshold;
P.reset_thr = -40; % Threshold to reset spike detection;
P.V0 = -63.6; % voltage initial condition in mV (chosen near Meng VU fixed point)
P.U0 = 0.43; % adaptation variable initial condition (chosen near Meng VU fixed point)


NS = runNoiseSims(P);

plotFirstSpike_boo = 1; % 1 to only plot the first spike per sim in the histograms, 
% 0 to plot all spikes recorded during each simulation.
if plotFirstSpike_boo
    spikes_toplot = NS.first_spikes';
else
    spikes_toplot = reshape([NS.spikes_cell{:}], [], 1);
end

% % Plot a single sim results
% figure;
% plot(P.tarray, NS.V_mat(:,1))
% hold on;
% plot(P.tarray, NS.U_mat(:,1).*20)
% plot(P.tarray, NS.input_mat(:,1)./100)
% legend(["V","20*U","I/100"])
% title(sprintf("Single Noisy Simulation \n A = %d, beta = %d", P.A, P.beta))
% xlabel("Time (ms)")

% Plot spike histogram
figure;
hold on;
% Histogram uses a linspace that spans from the start to the end of the
% input.
switch P.input_type
    case 'tent'
        h = histogram(spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(length(spikes_toplot)))));
        title(sprintf("%s Spiking Across %d Sims \n A = %d, beta = %d, sigma = %d",P.model, P.numSims, P.A, P.beta, P.sigma))
    case 'rectsin'
        h = histogram(spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(length(spikes_toplot)))));
        title(sprintf("%s Spiking Across %d Sims \n width = %d, beta = %d, sigma = %d",P.model, P.numSims, P.sin_width, P.beta, P.sigma))
end
plot(P.tarray, NS.input_noNoise.*(0.5*max(h.BinCounts)/max(NS.input_noNoise)))
legend(["Spike Counts", "I scaled"])
xlabel("Time (ms)")
ylabel("Number of Spikes")




%% Using dynamical threshold function

% Input params
P.t0 = 1; % start time of input. Make sure this is long enough to allow V and U 
% to reach steady states before input arrives.
P.inputShift = 5; % The shift in time between the start of the first input and 
% the start of the second input in a single period, when there are two
% inputs being provided. In ms.
P.isSecondInput = 1; % boolean for second input current profile.
P.A = 400; % Slope of the tent input.  
P.beta = 400; % Height of the tent input peak.
P.sin_width = 5;

% Simulation params struct
P.model = 'VU phasic';
P.input_type = 'rectsin';

switch P.input_type
    case 'tent'
        % Input as a tent function
        tent = @(t,t0,A,beta) (t<(t0+beta/A))*A*(t-t0) + (t>=(t0+beta/A))*(-A*(t-(t0+2*beta/A)));
        % The tent function drops down to negative instead of stopping at zero, so
        % we always need to take the maximum value between the tent function and 0
        % to prevent negative inputs.
        P.input = @(t,t0) (t>t0)*(t<=(t0+2*P.beta/P.A))*tent(t,t0,P.A,P.beta);
        P.Tlength = 2*P.beta/P.A + P.isSecondInput*(P.inputShift + 2*P.beta/P.A); % Length of a single period (ms).
    case 'rectsin'
        % input as a rectified sine wave. Width represents the width of the
        % single positive half of the wave, and so is half the length of
        % what the total period would be.
        sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
        P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
        P.Tlength = P.sin_width + P.isSecondInput*(P.inputShift + P.sin_width); % Length of a single period (ms).
    otherwise
        error("Invalid input type provided.")
end

P.dt = 0.001;
P.numT = 1; % Number of input periods to run per simulation.
P.simLength = P.t0 + P.numT*P.Tlength;
% Period start times. Each period begins a new set of input(s).
P.T_startTimes = P.t0 : P.Tlength : P.simLength - P.Tlength;
P.tarray = [0 : P.dt : P.simLength]';
P.numSteps = length(P.tarray);
P.spike_thr = 0; % Spike detection voltage threshold;
P.reset_thr = -40; % Threshold to reset spike detection;
P.V0 = -63.6; % voltage initial condition in mV (chosen near Meng VU fixed point)
P.U0 = 0.43; % adaptation variable initial condition (chosen near Meng VU fixed point)

% Kicking params
P.kickDensity = 10; % How many kicks to have per milisecond.
% KickIndices is a better replacement for kickTimes as the basis of
% determing when to kick, because converting from kickTimes to indices has
% rounding errors, but not when converting from indices to times. The
% density of kicks is P.kickDensity per milisecond. Note I don't
% kick all the way until the end of the simulation because at the very end
% it requires large kick amplitudes in order to get spikes to occur before
% the simulation ends, so it's best to just leave that out.
P.kickIndices = 1 : round((1/P.kickDensity)*(1/P.dt)) : P.numSteps - 0.5/P.dt;
% Evenly spaced kick times of P.kickDensity density, calculated from
% P.kickIndices.
P.kickTimes = P.kickIndices.*P.dt;
P.kickIncrement = 0.1; % Size to increase kick amplitude each iteration.


checkSimSpike(P); % checks that baseline simulation doesn't spike.
DT = getDynThr(P); % run dynamical threshold simulations.


% Get spike histogram to plot with DT
P.sigma = 100;
P.numSims = 10000; % Number of simulations to run.
NS = runNoiseSims(P); % run simulations with noise to get spiking distribution
% to plot alongside DT
plotFirstSpike_boo = 1; % 1 to only plot the first spike per sim in the histograms, 
% 0 to plot all spikes recorded during each simulation.
if plotFirstSpike_boo
    spikes_toplot = NS.first_spikes';
else
    spikes_toplot = reshape([NS.spikes_cell{:}], [], 1);
end

% and plot with DT
figure;
hold on;

switch P.input_type
    case 'tent'
        h = histogram(spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(2*length(spikes_toplot)))));
        title(sprintf("%s Dynamical Threshold \n A = %d, beta = %d",P.model, P.A, P.beta))
    case 'rectsin'
        h = histogram(spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(2*length(spikes_toplot)))), "FaceColor","y");
        title(sprintf("%s Dynamical Threshold \n width = %d, beta = %d",P.model, P.sin_width, P.beta))
end
plot(P.tarray, DT.input.*0.25*(max(h.BinCounts)/max(DT.input)),"color",'b')
plot(P.kickTimes, DT.kickSizes.*(max(h.BinCounts)/max(DT.kickSizes)), "color",'r')
plot(P.kickTimes, DT.spike_lags.*(max(h.BinCounts)/max(DT.spike_lags)), "color",'g')
xlabel("Time (ms)")
legend(["spike counts","I scaled","Dyn Thr scaled","spike lag scaled"])



% DT with adjusted spiketime histogram
% To account for the different spike time lags that occur at different
% points in the dynamical threshold, I need to subtract the proper spike
% lag from each spike in the distribution generated with noisy input. To do
% this, I need to interpolate the spike lag curve with the actual
% distribution spike times.

% spike lags for all spikes in the distribution
dist_spikelags = interp1(P.kickTimes, DT.spike_lags, spikes_toplot, 'linear', 'extrap');
% Subtract the time lags from each spike time
shifted_spikes_toplot = spikes_toplot - dist_spikelags;

figure;
hold on;
plot(P.kickTimes, DT.spike_lags,'color','g')
scatter(spikes_toplot,dist_spikelags,'MarkerEdgeColor','yellow')
legend(["Spike lag curve","Interpolation points"])
title("Checking Interpolation of Noisy Spiketime Distribution")
xlabel("kick and spike times (ms)")
ylabel("spike lag time (ms)")

figure;
hold on;

switch P.input_type
    case 'tent'
        h = histogram(shifted_spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(2*length(spikes_toplot)))), "FaceColor","magenta");
        title(sprintf("Dynamical Threshold \n A = %d, beta = %d", P.A, P.beta))
    case 'rectsin'
        h = histogram(shifted_spikes_toplot,...
            linspace(P.t0, P.simLength, 25+ceil(sqrt(2*length(spikes_toplot)))), "FaceColor","magenta");
        title(sprintf("Dynamical Threshold, Shifted Spike Histogram \n width = %d, beta = %d", P.sin_width, P.beta))
end
plot(P.tarray, DT.input.*0.25*(max(h.BinCounts)/max(DT.input)),"color",'b')
plot(P.kickTimes, DT.kickSizes.*(max(h.BinCounts)/max(DT.kickSizes)), "color",'r')
plot(P.kickTimes, DT.spike_lags.*(max(h.BinCounts)/max(DT.spike_lags)), "color",'g')
xlabel("Time (ms)") 
legend(["shifted spike counts","I scaled","Dyn Thr scaled","spike lag scaled"])



%% Phase plane analysis
% Run a single simulation and plot the results in the VU phase plane

% Input params
P.t0 = 1; % start time of input. Make sure this is long enough to allow V and U 
% to reach steady states before input arrives.
P.inputShift = 5; % The shift in time between the start of the first input and 
% the start of the second input in a single period, when there are two
% inputs being provided. In ms.
P.isSecondInput = 0; % boolean for second input current.
P.A = 800; % Slope of the tent input.  
P.beta = 400; % Height of the tent input peak.
P.sigma = 0; % Noise amplitude
P.sin_width = 1;

% Simulation params struct
P.model = 'VU tonic';
P.input_type = 'rectsin';
P.numSims = 1; % Number of simulations to run.
P.numT = 3; % Number of input periods to run per simulation.
P.dt = 0.001;
switch P.input_type
    case 'tent'
        % Input as a tent function
        tent = @(t,t0,A,beta) (t<(t0+beta/A))*A*(t-t0) + (t>=(t0+beta/A))*( -A*(t-(t0+2*beta/A)));
        % The tent function drops down to negative instead of stopping at zero, so
        % we always need to take the maximum value between the tent function and 0
        % to prevent negative inputs.
        P.input = @(t,t0) (t>t0)*(t<=(t0+2*P.beta/P.A))*tent(t,t0,P.A,P.beta);
        P.simLength = P.t0 + P.numT*4*P.beta/P.A;
        P.T_startTimes = P.t0 : 4*P.beta/P.A : P.simLength - 4*P.beta/P.A; % All period
        % start times.
    case 'rectsin'
        % input as a rectified sine wave. Width represents the width of the
        % single positive half of the wave, and so is half the length of
        % what the total period would be.
        sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
        P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
        P.simLength = P.t0 + P.numT*2*P.sin_width;
        P.T_startTimes = P.t0 : 2*P.sin_width : P.simLength - 2*P.sin_width;
    otherwise
        error("Invalid input type provided.")
end
P.tarray = [0 : P.dt : P.simLength]';
P.numSteps = length(P.tarray);
P.spike_thr = 0; % Spike detection voltage threshold;
P.reset_thr = -40; % Threshold to reset spike detection;
P.V0 = -63.6; % voltage initial condition in mV (chosen near Meng VU fixed point)
P.U0 = 0.43; % adaptation variable initial condition (chosen near Meng VU fixed point)


PPSim = runNoiseSims(P);

% Plot a single sim results
figure;
plot(P.tarray, PPSim.V_mat(:,1))
hold on;
plot(P.tarray, PPSim.U_mat(:,1).*20)
plot(P.tarray, PPSim.input_mat(:,1)./100)
legend(["V","20*U","I/100"])
title(sprintf("Single Simulation %s \n A = %d, beta = %d, sigma = %d",P.model, P.A, P.beta, P.sigma))
xlabel("Time (ms)")

% % Animate phase plane trajectory (slow)
% figure;
% title("VU Phase Plane Trajectory")
% xlabel("V (mV)")
% ylabel("U")
% for i = 1:P.numSteps
%     plot(PPSim.V_mat(1:i), PPSim.U_mat(1:i))
%     drawnow
% end

% Make line segments from the VU data for color plotting
Vseg = [PPSim.V_mat(1:end-1), PPSim.V_mat(2:end)];
Useg = [PPSim.U_mat(1:end-1), PPSim.U_mat(2:end)];

% Plot the segments just to see them
figure;
plt = plot(Vseg',Useg','-','LineWidth',4,'Visible','Off'); 
% axis equal;
xlim([min(PPSim.V_mat) max(PPSim.V_mat)]);
ylim([min(PPSim.U_mat) max(PPSim.U_mat)]);
switch P.input_type
    case 'tent'
        title(sprintf("Phase Plane Trajectory %s \n numT = %d, A = %d, beta = %d, sigma = %d",P.model, P.numT, P.A, P.beta, P.sigma))
    case 'rectsin'
        title(sprintf("Phase Plane Trajectory %s \n numT = %d, sinwidth = %d, beta = %d, sigma = %d",P.model, P.numT, P.sin_width, P.beta, P.sigma))
end
xlabel("V (mV)")
ylabel("U")

% Set all line segment colors
segColors = jet(size(Vseg,1)); % Choose a colormap
set(plt, {'Color'}, mat2cell(segColors,ones(size(Vseg,1),1),3))
set(plt, "Visible", "on");


%% Sigma vs time Spikes Heatmap
% I am going to create a 2D heatmap showing the distribution of spikes as
% functions of time (x-axis) and noise amplitude sigma (y-axis). 

% Input params
P.t0 = 1; 
P.inputShift = 5; 
P.isSecondInput = 1; 
P.A = 400; 
P.beta = 400; 
P.sin_width = 5;

% Simulation params struct
P.model = 'VU phasic';
P.input_type = 'rectsin';

switch P.input_type
    case 'tent'
        % Input as a tent function
        tent = @(t,t0,A,beta) (t<(t0+beta/A))*A*(t-t0) + (t>=(t0+beta/A))*(-A*(t-(t0+2*beta/A)));
        P.input = @(t,t0) (t>t0)*(t<=(t0+2*P.beta/P.A))*tent(t,t0,P.A,P.beta);
        P.Tlength = 2*P.beta/P.A + P.isSecondInput*(P.inputShift + 2*P.beta/P.A); % Length of a single period (ms).
    case 'rectsin'
        sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
        P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
        P.Tlength = P.sin_width + P.isSecondInput*(P.inputShift + P.sin_width); % Length of a single period (ms).
    otherwise
        error("Invalid input type provided.")
end

P.numSims = 1000; 
P.dt = 0.001;
P.numT = 1;
P.simLength = P.t0 + P.numT*P.Tlength;
P.T_startTimes = P.t0 : P.Tlength : P.simLength - P.Tlength;
P.tarray = [0 : P.dt : P.simLength]';
P.numSteps = length(P.tarray);
P.spike_thr = 0; 
P.reset_thr = -40; 
P.V0 = -63.6; 
P.U0 = 0.43;

hist_edges = linspace(P.t0, P.t0 + P.Tlength, 25+ceil(sqrt(P.numSims)));

P.sigma = 0;
checkSimSpike(P); % checks that baseline simulation doesn't spike.


sigma_range = 0 : 5 : 250; % Define sigma values to use.

counts_mat = zeros(length(sigma_range), length(hist_edges)-1); % Matrix to hold spike count data.


for i = 1:length(sigma_range)
    
    P.sigma = sigma_range(i); % Assign sigma.

    NS = runNoiseSims(P); % Get the data.
    
    plotFirstSpike_boo = 1; % 1 to only plot/use the first spike per sim in the histograms, 
    % 0 to plot/use all spikes recorded during each simulation.
    if plotFirstSpike_boo
        spikes_toplot = NS.first_spikes';
    else
        spikes_toplot = reshape([NS.spikes_cell{:}], [], 1);
    end

    [counts, edges] = histcounts(spikes_toplot, hist_edges); % Bin the data.
    counts_mat(i,:) = counts;
    
    fprintf("sigma parameter %d out of %d completed \n", i, length(sigma_range))
end



%% plot
figure;
hold on;
[X,Y] = meshgrid(hist_edges(2:end),sigma_range);
im = imagesc(hist_edges(2:end),sigma_range,counts_mat);
cb = colorbar;
cb.Label.String = "Number of Spikes";
switch P.input_type
    case 'tent'
        title(sprintf("Spike Count Image tent input \n A = %d, beta = %d, IS = %d", P.A, P.beta, P.inputShift))
    case 'rectsin'
        title(sprintf("Spike Count Image rectsin input \n width = %d, beta = %d, IS = %d", P.sin_width, P.beta, P.inputShift))
end
xlabel("Time (ms)")
ylabel("Sigma")
% axis image;
% ax = gca;
% ax.YTickLabel = linspace(min(sigma_range), max(sigma_range), numel(ax.YTick)+1);
% ax.XTickLabel = round(linspace(min(hist_edges), max(hist_edges), numel(ax.XTick)+1), 2);

plot(P.tarray, NS.input_noNoise.*(size(counts_mat,1)/max(NS.input_noNoise)), "color", 'k')

% 
% figure;
% hold on;
% plot(sigma_range, sum(counts_mat,2))
% xlabel("sigma")
% ylabel("Total Spike Count")
% title("Total Number of Spikes in Histograms")

%% Plot same heatmap as above but with dynamical threshold trace included for the noiseless case

% Kicking params
P.kickDensity = 5; % How many kicks to have per milisecond.
% Evenly spaced kick times of P.kickDensity density. Note kickTimes can't
% start at 0 because I get an indexing by 0 error later. Also note I don't
% kick all the way until the end of the simulation because at the very end
% it requires large kick amplitudes in order to get spikes to occur before
% the simulation ends.
P.kickTimes = linspace(P.dt, (P.simLength - 0.5), P.kickDensity*((P.simLength - 0.5) - P.dt));
P.kickIncrement = 0.1; % Size to increase kick amplitude each iteration.

checkSimSpike(P); % checks that baseline simulation doesn't spike.
DT = getDynThr(P); % run dynamical threshold simulations.

% Plot
figure;
hold on;
im = imagesc(counts_mat);
cb = colorbar;
cb.Label.String = "Number of Spikes";
switch P.input_type
    case 'tent'
        title(sprintf("Spike Count Image tent input \n A = %d, beta = %d, IS = %d", P.A, P.beta, P.inputShift))
    case 'rectsin'
        title(sprintf("Spike Count Image rectsin input \n width = %d, beta = %d, IS = %d", P.sin_width, P.beta, P.inputShift))
end
xlabel("Time (ms)")
ylabel("Sigma")
axis image;
ax = gca;
ax.YTickLabel = linspace(min(sigma_range), max(sigma_range), numel(ax.YTick)+1);
ax.XTickLabel = round(linspace(0, max(hist_edges), numel(ax.XTick)+1), 2);

plot(P.tarray.*(size(counts_mat,2)/max(P.tarray)), NS.input_noNoise.*(size(counts_mat,1)/max(NS.input_noNoise)), "color", 'k')
plot(P.kickTimes.*(size(counts_mat,2)/max(P.kickTimes)), DT.kickSizes.*(size(counts_mat,1)/max(DT.kickSizes)), "color",'r') % plot dynamical threshold

%% 3D heatmap with input amplitude on the third dimension **Note: Only works for tent currently**
% This will be an extension of the 2D heatmap above. On the third dimension
% the amplitude of the input current provided will be varied. There are
% three possible ways to do this: 1) increasing slope and keeping width the
% same 2) increasing width and keeping slope the same 3) increasing height 
% and keeping area the same. For the sake of simplicity I will start by
% keeping width the same and increasing height/slope, that way all the
% simulation times are the same. 

% Input params
P.t0 = 1; % start time of input. Make sure this is long enough to allow V and U 
% to reach steady states before input arrives.
P.inputShift = 1; % The shift in time between the start of the first input and 
% the start of the second input in a single period, when there are two
% inputs being provided. In ms.
P.isSecondInput = 0; % boolean for second input current.

% Simulation params struct
P.model = 'VU phasic';
P.input_type = 'tent';
P.numSims = 200; % Number of simulations to run per parameter set.
P.dt = 0.001;
P.spike_thr = 0; % Spike detection voltage threshold;
P.reset_thr = -40; % Threshold to reset spike detection;
P.V0 = -63.6; % voltage initial condition in mV (chosen near Meng VU fixed point)
P.U0 = 0.43; % adaptation variable initial condition (chosen near Meng VU fixed point)


% Define a set width for the input current profiles (in ms).
input_width = 4;

beta_range = 900 : 25 : 1300; % input height values to use.

sigma_range = 0 : 0.2 : 15; % Define sigma values to use.

% The binning to use for the heatmap. Needs to be the same for all sigma
% values.
hist_edges = linspace(P.t0, P.t0 + input_width, 25+ceil(sqrt(P.numSims)));

counts_mat = zeros(length(sigma_range), length(hist_edges)-1, length(beta_range)); % Matrix to hold spike count data.


for j = 1:length(beta_range)

    P.beta = beta_range(j);
    P.A = 2*P.beta/input_width; % This should keep the total width of the input
    % the same across all beta values.

    switch P.input_type
        case 'tent'
            % Input as a tent function
            tent = @(t,t0,A,beta) beta*(t-t0)*(t>t0)-2*beta*(t-t0-A/beta)*(t>(t0+A/beta));
            % The tent function drops down to negative instead of stopping at zero, so
            % we always need to take the maximum value between the tent function and 0
            % to prevent negative inputs.
            P.input = @(t,t0) max(tent(t,t0,P.A,P.beta),0);
            P.simLength = P.t0 + 4*P.beta/P.A;
        case 'rectsin'
            % input as a rectified sine wave. Width represents the width of the
            % single positive half of the wave, and so is half the length of
            % what the total period would be.
            sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
            P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
            P.simLength = P.t0 + 2*P.sin_width;
        otherwise
            error("Invalid input type provided.")
    end
    P.tarray = [0 : P.dt : P.simLength]';
    P.numSteps = length(P.tarray);

%     checkSimSpike(P); % Checks that the baseline simulation doesn't cause spiking.
    
    
    for i = 1:length(sigma_range)
        
        P.sigma = sigma_range(i); % Assign sigma.
    
        NS = runNoiseSims(P); % Get the data.
    
        [counts, edges] = histcounts(vertcat(NS.spikes_cell{:}), hist_edges); % Bin the data.
        counts_mat(i,:,j) = counts;
    
    end

    fprintf("beta parameter %d out of %d completed \n", j, length(beta_range))

end


%Plot 3D pseudocolor plot
[X,Y,Z] = meshgrid(hist_edges(1:end-1), sigma_range, beta_range);
pc3 = pcolor3(X,Y,Z,counts_mat, 'alpha', 0.1,'edgealpha', 0.05, 'alphalim',[0,P.numSims]);
xlabel("Time (ms)")
ylabel("Sigma")
zlabel("beta")
title("Beta vs Sigma Spike Histograms")
cb = colorbar;
cb.Label.String = "Number of Spikes";



%% Testing Different Spike Detection Thresholds
% The exact spike detection threshold may be impacting significantly the
% time at which spikes are recorded, possibly accounting for the strange
% fact that we see a higher probability of spiking after the minimum of the
% DT. 
% To test this, I am going to record the spike times of a noiseless
% simulation and observe how that changes as the detection threshold
% P.spike_thr changes. Then, I am going to add noise in and look at how the
% distribution of spikes around the DT valley changes as the threshold
% changes, if at all.

% Input params
P.t0 = 1; 
P.inputShift = 2; 
P.isSecondInput = 1; 
P.A = 400; 
P.beta = 1000; 
P.sin_width = 5;
P.sigma = 0;

% Simulation params struct
P.model = 'VU phasic';
P.input_type = 'rectsin';

switch P.input_type
    case 'tent'
        % Input as a tent function
        tent = @(t,t0,A,beta) (t<(t0+beta/A))*A*(t-t0) + (t>=(t0+beta/A))*(-A*(t-(t0+2*beta/A)));
        P.input = @(t,t0) (t>t0)*(t<=(t0+2*P.beta/P.A))*tent(t,t0,P.A,P.beta);
        P.Tlength = 2*P.beta/P.A + P.isSecondInput*(P.inputShift + 2*P.beta/P.A); % Length of a single period (ms).
    case 'rectsin'
        sinf = @(t,t0,beta,width) beta*sin((t-t0)*pi/width);
        P.input = @(t,t0) (t>t0)*(t<(t0+P.sin_width))*sinf(t,t0,P.beta,P.sin_width);
        P.Tlength = P.sin_width + P.isSecondInput*(P.inputShift + P.sin_width); % Length of a single period (ms).
    otherwise
        error("Invalid input type provided.")
end

P.numSims = 1; 
P.dt = 0.001;
P.numT = 1;
P.simLength = P.t0 + P.numT*P.Tlength;
P.T_startTimes = P.t0 : P.Tlength : P.simLength - P.Tlength;
P.tarray = [0 : P.dt : P.simLength]';
P.numSteps = length(P.tarray);
P.reset_thr = -40; 
P.V0 = -63.6; 
P.U0 = 0.43;

thr_range = -20 : 1 : 10; % the threshold values to use.
det_spike_times = zeros(size(thr_range)); % deterministic spike times for the
% no noise simulations.

% Note P.numSims should be 1 and P.sigma should be 0.
if P.numSims ~= 1 || P.sigma ~= 0
    error("P.numSims must be 1 and P.sigma must be 0 for this analysis.")
end

for i = 1:length(thr_range)
    P.spike_thr = thr_range(i);
    NS = runNoiseSims(P);

    if isempty(reshape([NS.spikes_cell{:}], [], 1))
        error("No spikes recorded. Change input to induce spikes without noise.")
    end

    det_spike_times(i) = reshape([NS.spikes_cell{:}], [], 1);
end

SDf = figure;
hold on;
plot(thr_range, det_spike_times)
title("Spike Detection Times")
xlabel("Threshold Value (mV)")
ylabel("Time of spike (ms)")


% Now adding noise and running many simulations to get an average for the
% spike times

P.sigma = 100;
P.numSims = 100;
thr_range = -20 : 0.1 : 10; % the threshold values to use.
noise_spike_times = zeros(size(thr_range));

for i = 1:length(thr_range)
    P.spike_thr = thr_range(i);
    NS = runNoiseSims(P);

    spikes_cell = NS.spikes_cell; % Extracts spike times from all simulations.
    % Some simulations have no spikes, and some have more than one.

    nonempty_inds = ~cellfun('isempty',spikes_cell); % Identify indices of sims 
    % that have a spike.
    first_spikes = cellfun(@(c) c(1), spikes_cell(nonempty_inds)); % Extract 
    % only the first spike time from all cells.

    noise_spike_times(i) = mean(first_spikes);

    fprintf("Threshold value %d out of %d completed.", i,length(thr_range))
end

set(0, 'currentfigure', SDf);
hold on;
plot(thr_range, noise_spike_times)
% title("Noisy Spike Detection Times")
% xlabel("Threshold Value (mV)")
% ylabel("Time of spike (ms)")
legend(["sigma=0", sprintf("sigma=%d",P.sigma)])
