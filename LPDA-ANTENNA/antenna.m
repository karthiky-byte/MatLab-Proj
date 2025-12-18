%% FINAL SMOOTHED DESIGN: 8-Element LPDA
%  Goal: Minimize Ripple, Maximize Gain, Show Dimensions & Performance
clear; clc; close all;

%% === 1. OPTIMIZED PARAMETERS ===
Tune_Boom_Limit = 0.60;   
Tune_Sigma      = 0.14;   
Tune_Tau        = 0.87;   
Tune_Wire_Width = 0.005;  
f_min = 500e6;
f_max = 2e9;
f_center = 1.25e9; 
c = 3e8;

%% === 2. GEOMETRY GENERATION ===
L_max = (c / f_min) / 2; 
N = 8; 
% Generate Lengths
L = zeros(1, N); L(1) = L_max;
for i = 2:N, L(i) = L(i-1) * Tune_Tau; end

% Generate Spacing
S = zeros(1, N-1);
for i = 1:N-1, S(i) = 2 * Tune_Sigma * L(i+1); end

% Scale to fit 600mm
if sum(S) > Tune_Boom_Limit
    scale = Tune_Boom_Limit / sum(S) * 0.99;
    S = S * scale;
end
W = ones(1, N) * Tune_Wire_Width;

%% === 3. BUILD ANTENNA ===
ant = lpda;
ant.ArmLength = L;
ant.ArmSpacing = S;      
ant.ArmWidth = W;
ant.FeedWidth = 0.004; 

% Air Substrate
try ant.Substrate = dielectric('Name', 'Air', 'EpsilonR', 1, 'LossTangent', 0); catch; ant.Substrate = []; end
try ant.Length = 1.2; ant.Width = 1.2; catch; try ant.BoardLength = 1.2; ant.BoardWidth = 1.2; catch; end; end

% Mesh (Silent)
mesh(ant, 'MaxEdgeLength', (c/f_max)/7); 

%% === 4. S11 SIMULATION ===
fprintf('Simulating S11 and Gain... Please wait.\n');
freq_range = linspace(f_min, f_max, 200); 
S_param = sparameters(ant, freq_range);
S11dB = 20*log10(abs(squeeze(S_param.Parameters(1,1,:))));

% Plot S11
figure('Name', 'Return Loss');
plot(freq_range/1e9, S11dB, 'b-', 'LineWidth', 2);
yline(-10, 'r--', 'Limit (-10dB)');
grid on; xlabel('Frequency (GHz)'); ylabel('S11 (dB)');
title('Return Loss (S11)');

%% === 5. GAIN CALCULATION ===
GainCurve = zeros(size(freq_range));
for k = 1:length(freq_range)
    GainCurve(k) = max(max(pattern(ant, freq_range(k)))); 
end

% Plot Gain
figure('Name', 'Gain vs Frequency');
plot(freq_range/1e9, GainCurve, 'r-', 'LineWidth', 2);
yline(6, 'k--', 'Target 6dBi');
grid on; xlabel('Frequency (GHz)'); ylabel('Gain (dBi)');
title('Realized Gain');

%% === 6. RESULT ANALYSIS & PRINTING ===
% Find Min S11 (Best Matching)
[min_s11, idx_s11] = min(S11dB);
freq_s11 = freq_range(idx_s11);

% Find Max Gain (Best Radiation)
[max_gain, idx_gain] = max(GainCurve);
freq_gain = freq_range(idx_gain);

clc; % Clear screen to show only the final report

% --- PRINT GEOMETRY TABLE ---
fprintf('\n======================================================\n');
fprintf('           ANTENNA GEOMETRY (Dimensions)              \n');
fprintf('======================================================\n');
fprintf(' El # | Length (mm) | Spacing to Next (mm) | Width (mm)\n');
fprintf('------------------------------------------------------\n');
for i = 1:N
    if i < N
        spacing_val = S(i) * 1000; % Convert to mm
    else
        spacing_val = 0;           % Last element has no next spacing
    end
    fprintf('  %2d  |   %6.2f    |       %6.2f         |   %5.2f\n', ...
        i, L(i)*1000, spacing_val, W(i)*1000);
end
fprintf('------------------------------------------------------\n');
fprintf(' TOTAL BOOM LENGTH: %.2f mm\n', sum(S)*1000);
fprintf('======================================================\n\n');

% --- PRINT PERFORMANCE SUMMARY ---
fprintf('============================================\n');
fprintf('          FINAL SIMULATION RESULTS          \n');
fprintf('============================================\n');
fprintf('FREQUENCY RANGE: %.1f MHz - %.1f GHz\n', f_min/1e6, f_max/1e9);
fprintf('--------------------------------------------\n');
fprintf('MAX GAIN       : %.2f dBi  @ %.2f GHz\n', max_gain, freq_gain/1e9);
fprintf('MIN RETURN LOSS: %.2f dB   @ %.2f GHz\n', min_s11, freq_s11/1e9);
fprintf('============================================\n\n');

%% === 7. 3D RADIATION PATTERN ===
figure('Name','3D Radiation Pattern');
pattern(ant, f_center);
title(['3D Radiation Pattern @ ' num2str(f_center/1e9) ' GHz']);