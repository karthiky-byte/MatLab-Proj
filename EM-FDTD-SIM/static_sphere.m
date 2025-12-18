%% Realistic Doppler Effect in 2D FDTD with Moving Object
clear; clc;

% ---------------- Grid & Physical Parameters ----------------
Nx = 200; Ny = 200;
dx = 1e-3; dy = 1e-3;
c0 = 3e8; eps0 = 8.854e-12; mu0 = 4*pi*1e-7;
dt = 1/(c0*sqrt(1/dx^2 + 1/dy^2)) * 0.99;
n_steps = 1000;

% Initialize fields
Ez = zeros(Nx, Ny);
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);

% Relative permittivity & permeability
eps_r = ones(Nx, Ny);
mu_r = ones(Nx, Ny);

% PML (Simple absorbing boundary)
pml_size = 10;
sigma_max = 1e3;
sigma = zeros(Nx, Ny);
sigma(1:pml_size,:) = sigma_max;
sigma(end-pml_size+1:end,:) = sigma_max;
sigma(:,1:pml_size) = sigma_max;
sigma(:,end-pml_size+1:end) = sigma_max;

% ---------------- Source Parameters ----------------
source_x = 50; source_y = 100;
f0 = 1e9;      % source frequency (Hz)
t0 = 50; spread = 20;  % Gaussian envelope

% ---------------- Moving Dielectric Object ----------------
moving_objects = [
    struct('x',120,'y',100,'w',12,'h',8,'dx',0.02,'dy',0.0,'eps',6.0)
];

% ---------------- Doppler Effect Probes ----------------
probe_front = [150,100];  % in front of moving object
probe_back  = [30,100];   % behind moving object
Ez_front = zeros(1, n_steps);
Ez_back  = zeros(1, n_steps);

% ---------------- Main FDTD Loop ----------------
figure;
for n = 1:n_steps
    % ---- Reset permittivity ----
    eps_r(:) = 1;

    % ---- Update moving objects ----
    for i = 1:length(moving_objects)
        obj = moving_objects(i);
        x0 = max(1, round(obj.x)); y0 = max(1, round(obj.y));
        x1 = min(round(obj.x + obj.w - 1), Nx);
        y1 = min(round(obj.y + obj.h - 1), Ny);
        eps_r(x0:x1, y0:y1) = obj.eps;

        % Update position & bounce
        moving_objects(i).x = obj.x + obj.dx;
        moving_objects(i).y = obj.y + obj.dy;
        if moving_objects(i).x+obj.w>=Nx-pml_size || moving_objects(i).x<=1+pml_size
            moving_objects(i).dx = -moving_objects(i).dx;
        end
    end

    % ---- Update Magnetic Fields ----
    Hx(:,1:end-1) = (1 - dt*sigma(:,1:end-1)./(2*mu0)) ./ (1 + dt*sigma(:,1:end-1)./(2*mu0)) .* Hx(:,1:end-1) ...
        - (dt ./ (mu0 * mu_r(:,1:end-1) * dy)) ./ (1 + dt*sigma(:,1:end-1)./(2*mu0)) .* diff(Ez,1,2);

    Hy(1:end-1,:) = (1 - dt*sigma(1:end-1,:)./(2*mu0)) ./ (1 + dt*sigma(1:end-1,:)./(2*mu0)) .* Hy(1:end-1,:) ...
        + (dt ./ (mu0 * mu_r(1:end-1,:) * dx)) ./ (1 + dt*sigma(1:end-1,:)./(2*mu0)) .* diff(Ez,1,1);

    % ---- Curl of H ----
    curlH = (Hy(2:end-1,2:end-1) - Hy(1:end-2,2:end-1)) / dx ...
          - (Hx(2:end-1,2:end-1) - Hx(2:end-1,1:end-2)) / dy;

    % ---- Update Electric Field ----
    c1 = (1 - dt*sigma(2:end-1,2:end-1)./(2*eps0)) ./ (1 + dt*sigma(2:end-1,2:end-1)./(2*eps0));
    c2 = (dt ./ (eps0 * eps_r(2:end-1,2:end-1))) ./ (1 + dt*sigma(2:end-1,2:end-1)./(2*eps0));
    Ez(2:end-1,2:end-1) = c1 .* Ez(2:end-1,2:end-1) + c2 .* curlH;

    % ---- Source Excitation ----
    Ez(source_x, source_y) = Ez(source_x, source_y) + ...
        2.0 * exp(-((n-t0)/spread)^2) .* cos(2*pi*f0*n*dt);

    % ---- Record Probe Data (for Doppler Effect) ----
    Ez_front(n) = Ez(probe_front(1), probe_front(2));
    Ez_back(n)  = Ez(probe_back(1), probe_back(2));

    % ---- Print intensity every 50 steps ----
    if mod(n,50)==0
        H_intensity = sqrt(Hx.^2 + Hy.^2);
        fprintf('Step %d: Ez[max]=%.3e |H|max=%.3e\n', ...
            n, max(abs(Ez(:))), max(H_intensity(:)));
    end

    % ---- Visualization (Ez & |H|) ----
    if mod(n,20)==0
        H_intensity = sqrt(Hx.^2 + Hy.^2);

        subplot(1,2,1);
        imagesc(Ez'); axis equal tight;
        title(['Ez Field at Step ', num2str(n)]);
        xlabel('X'); ylabel('Y');
        colorbar; colormap jet; caxis([-0.05 0.05]);
        hold on;
        for i = 1:length(moving_objects)
            obj = moving_objects(i);
            rectangle('Position',[obj.x obj.y obj.w obj.h],'EdgeColor','white');
        end
        hold off;

        subplot(1,2,2);
        imagesc(H_intensity'); axis equal tight;
        title('|H| Magnetic Intensity');
        xlabel('X'); ylabel('Y');
        colorbar; colormap hot; caxis auto;
        hold on;
        for i = 1:length(moving_objects)
            obj = moving_objects(i);
            rectangle('Position',[obj.x obj.y obj.w obj.h],'EdgeColor','cyan');
        end
        hold off;

        drawnow;
    end
end

%% ---------------- Doppler Effect Analysis ----------------
fs = 1/dt;
N = n_steps;
f_axis = (0:N-1) * (fs/N);
Ez_front_fft = abs(fft(Ez_front));
Ez_back_fft  = abs(fft(Ez_back));

% Extract frequency near source f0
[~, idx_front] = max(Ez_front_fft(round(f0*N/fs*0.5):round(f0*N/fs*1.5)));
idx_front = idx_front + round(f0*N/fs*0.5) - 1;
[~, idx_back]  = max(Ez_back_fft(round(f0*N/fs*0.5):round(f0*N/fs*1.5)));
idx_back  = idx_back + round(f0*N/fs*0.5) - 1;

fprintf('\n--- Doppler Effect (Realistic) ---\n');
fprintf('Source f0 ≈ %.2f GHz\n', f0/1e9);
fprintf('Detected Front Frequency ≈ %.4f GHz\n', f_axis(idx_front)/1e9);
fprintf('Detected Back Frequency  ≈ %.4f GHz\n', f_axis(idx_back)/1e9);

% Plot Time-domain Signals
figure;
plot((1:n_steps)*dt, Ez_front,'r'); hold on;
plot((1:n_steps)*dt, Ez_back,'b');
legend('Front (Approaching)','Back (Receding)');
xlabel('Time (s)'); ylabel('Ez(t)');
title('Time-domain Electric Field at Probes');

% Plot Frequency-domain (Doppler Shift)
figure;
plot(f_axis(1:N/2), Ez_front_fft(1:N/2),'r','LineWidth',1.5); hold on;
plot(f_axis(1:N/2), Ez_back_fft(1:N/2),'b','LineWidth',1.5);
xlim([0 2*f0]);  % zoom around source frequency
legend('Front (Approaching)','Back (Receding)');
xlabel('Frequency (Hz)'); ylabel('|Ez(f)|');
title('Doppler Effect: Frequency Shift');
grid on;