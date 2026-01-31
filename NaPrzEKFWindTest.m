%% EKF MONTE CARLO: IMU + GPS(p) + GPS(v) + PITOT(q) + wind state
clc; clear; close all;

%% 1) SETTINGS
dt_imu = 0.01;
T      = 60;
t      = (0:dt_imu:T)';
N      = numel(t);
acc    = 0.6;

dt_gps   = 0.2;
gps_step = round(dt_gps/dt_imu);

dt_pitot   = 0.005;
pitot_step = max(1, round(dt_pitot/dt_imu));

rho = 1.225;

Nmc = 500;   % <-- set how many runs you want

%% 2) "TRUTH" (deterministic) will be the basis for the sensors to add their
%% noise and give to the filter.
a_true = zeros(N,1);
a_true(t>=2  & t<8)  =  1.2*acc;
a_true(t>=15 & t<18) = -1.5*acc;
a_true(t>=25 & t<35) =  0.6*acc;
a_true(t>=45 & t<50) = -0.8*acc;

% True wind: 0->2 m/s in 10s and back in 10s(very simplified model, but it
% doesnt really matter from what i have read).
tau    = mod(t, 20);
w_true = 1 - cos(pi * tau/10);

% Ground truth integration
v_true = zeros(N,1);
p_true = zeros(N,1);
for k = 2:N
    v_true(k) = v_true(k-1) + a_true(k-1)*dt_imu;
    p_true(k) = p_true(k-1) + v_true(k-1)*dt_imu + 0.5*a_true(k-1)*dt_imu^2;
end

v_air_true = v_true - w_true;
q_true     = 0.5*rho*(v_air_true.^2);

% Measurement time flags (fixed)
gps_new = false(N,1);   gps_new(1:gps_step:end) = true;
pitot_new = false(N,1); pitot_new(1:pitot_step:end) = true;

%% 3) NOISE / FILTER KNOBS
% IMU
sigma_a_meas = 0.20;
sigma_b_rw   = 0.01;

% GPS
sigma_gps_pos = 1.0;

% "GPS velocity" (pseudo Doppler, when tried discrete integration the gps
% noise destroyed the velocity readings)
sigma_gps_vel = 0.4;
R_gps_v = sigma_gps_vel^2;
H_gps_v = [0 1 0 0];

% PITOT q noise
dq_3sigma = 1.0;
sigma_q   = dq_3sigma/3;
R_pitot   = sigma_q^2;

% Process noise
q_a = (0.20)^2;
q_b = (sigma_b_rw^2)*dt_imu;

% Wind Gauss–Markov (filter model taken straight from internet)
tau_w = 10;                         % [s]
alpha = exp(-dt_imu/tau_w);
sigma_w_ss = 1.0;                   % [m/s]
q_w = (1 - alpha^2) * sigma_w_ss^2;

% Gate / pitot trust vs low airspeed
vmin_pitot = 1.5;
v0 = 2.0;

% GPS velocity smoothing
beta = 0.4;  % 0..1

%% 4) EKF MATRICES (constant)
F = [1 dt_imu -0.5*dt_imu^2  0;
     0 1     -dt_imu         0;
     0 0      1              0;
     0 0      0              alpha];

B_u = [0.5*dt_imu^2; dt_imu; 0; 0];

H_gps = [1 0 0 0];
R_gps = sigma_gps_pos^2;

Q = [0.25*dt_imu^4*q_a, 0.5*dt_imu^3*q_a, 0, 0;
     0.5*dt_imu^3*q_a,  dt_imu^2*q_a,     0, 0;
     0,                0,                q_b, 0;
     0,                0,                0,   q_w];

I = eye(4);

%% 5) STORAGE
MAE_p_all     = zeros(Nmc,1);
MAE_v_all     = zeros(Nmc,1);
MAE_w_all     = zeros(Nmc,1);
MAE_b_all     = zeros(Nmc,1);
MAE_vair_all  = zeros(Nmc,1);

% Optional: keep one “example run” trajectory to plot afterwards used
% mainly for deducing what went wrong in the filter
keepExample = true;
ex_idx = 1;
xhat_hist_ex = [];
p_gps_ex = [];
q_pitot_ex = [];

%% 6) PARFOR fallback
usePar = ~isempty(gcp('nocreate'));

if usePar
    parfor mc = 1:Nmc
        [MAE_p_all(mc), MAE_v_all(mc), MAE_w_all(mc), MAE_b_all(mc), MAE_vair_all(mc)] = ...
            run_one_EKF(dt_imu, N, a_true, p_true, v_true, w_true, q_true, ...
                        gps_new, pitot_new, ...
                        sigma_a_meas, sigma_b_rw, ...
                        sigma_gps_pos, ...
                        R_gps, H_gps, ...
                        beta, dt_gps, gps_step, R_gps_v, H_gps_v, ...
                        R_pitot, rho, vmin_pitot, v0, ...
                        F, B_u, Q, I);
    end
else
    for mc = 1:Nmc
        [MAE_p_all(mc), MAE_v_all(mc), MAE_w_all(mc), MAE_b_all(mc), MAE_vair_all(mc), ...
         xhat_hist_tmp, p_gps_tmp, q_pitot_tmp] = ...
            run_one_EKF(dt_imu, N, a_true, p_true, v_true, w_true, q_true, ...
                        gps_new, pitot_new, ...
                        sigma_a_meas, sigma_b_rw, ...
                        sigma_gps_pos, ...
                        R_gps, H_gps, ...
                        beta, dt_gps, gps_step, R_gps_v, H_gps_v, ...
                        R_pitot, rho, vmin_pitot, v0, ...
                        F, B_u, Q, I);

        if keepExample && mc == ex_idx
            xhat_hist_ex = xhat_hist_tmp;
            p_gps_ex = p_gps_tmp;
            q_pitot_ex = q_pitot_tmp;
        end
    end
end

%% 7) STATS
fprintf('\n MONTE CARLO STATS (%d runs) \n', Nmc);

print_stats("p [m]"     , MAE_p_all);
print_stats("v [m/s]"   , MAE_v_all);
print_stats("w [m/s]"   , MAE_w_all);
print_stats("b [m/s^2]" , MAE_b_all);
print_stats("v_air [m/s]", MAE_vair_all);

%% 8) PLOT ONE EXAMPLE RUN (only if ran with for-loop or saved elsewhere)
if keepExample && ~isempty(xhat_hist_ex)
    v_air_est_ex = xhat_hist_ex(:,2) - xhat_hist_ex(:,4);

    figure; grid on; hold on;
    plot(t, v_true, 'LineWidth', 1.2);
    plot(t, xhat_hist_ex(:,2), '--', 'LineWidth', 1.2);
    xlabel('t [s]'); ylabel('v [m/s]');
    legend('v true','v est'); title('Example run: v');

    figure; grid on; hold on;
    plot(t, w_true, 'LineWidth', 1.2);
    plot(t, xhat_hist_ex(:,4), '--', 'LineWidth', 1.2);
    xlabel('t [s]'); ylabel('w [m/s]');
    legend('w true','w est'); title('Example run: wind');

    figure; grid on; hold on;
    plot(t, p_true, 'LineWidth', 1.2);
    plot(t, xhat_hist_ex(:,1), '--', 'LineWidth', 1.2);
    plot(t, p_gps_ex, ':', 'LineWidth', 1.2);
    xlabel('t [s]'); ylabel('p [m]');
    legend('p true','p est','p gps'); title('Example run: position');

    figure; grid on; hold on;
    plot(t, q_true, 'LineWidth', 1.2);
    plot(t, q_pitot_ex, '.', 'MarkerSize', 5);
    xlabel('t [s]'); ylabel('q [Pa]');
    legend('q true','q pitot'); title('Example run: pitot q');

    figure; grid on; hold on;
    plot(t, v_air_true, 'LineWidth', 1.2);
    plot(t, v_air_est_ex, '--', 'LineWidth', 1.2);
    xlabel('t [s]'); ylabel('v_{air} [m/s]');
    legend('v_{air} true','v_{air} est'); title('Example run: v_air');
end

%% LOCAL FUNCTIONS
function print_stats(name, x)
    fprintf('\n%s\n', name);
    fprintf('  median(MAE): %.4f\n', median(x));
    fprintf('  IQR (25-75): %.4f .. %.4f\n', prctile(x,25), prctile(x,75));
    fprintf('  p95        : %.4f\n', prctile(x,95));
    fprintf('  std        : %.4f\n', std(x)); % sample std by default
end

function [MAE_p, MAE_v, MAE_w, MAE_b, MAE_vair, xhat_hist, p_gps, q_pitot] = run_one_EKF( ...
        dt_imu, N, a_true, p_true, v_true, w_true, q_true, ...
        gps_new, pitot_new, ...
        sigma_a_meas, sigma_b_rw, ...
        sigma_gps_pos, ...
        R_gps, H_gps, ...
        beta, dt_gps, gps_step, R_gps_v, H_gps_v, ...
        R_pitot, rho, vmin_pitot, v0, ...
        F, B_u, Q, I)

    % ---- IMU noise + bias RW ----
    b_true = zeros(N,1);
    for k = 2:N
        b_true(k) = b_true(k-1) + sigma_b_rw*sqrt(dt_imu)*randn;
    end
    a_imu = a_true + b_true + sigma_a_meas*randn(N,1);

    %  GPS position
    p_gps = p_true + sigma_gps_pos*randn(N,1);

    %  GPS pseudo velocity difference + low pass filter
    v_gps = nan(N,1);
    v_gps_f = v_gps;
    beta = 0.1; % 0..1 (Lower = Stronger smoothing)
    for k = 2:N
        if ~isnan(v_gps(k)) && ~isnan(v_gps_f(k-1))
            v_gps_f(k) = (1-beta)*v_gps_f(k-1) + beta*v_gps(k); % Doppler algorithm for smoother velocty calculation
        end
    end

    % ---- Pitot q ----
    q_pitot = q_true + sqrt(R_pitot)*randn(N,1);

    % ---- EKF init ----
    xhat = [p_gps(1); 0; 0; 0];
    P    = diag([25, 25, 0.01, 1.0^2]);

    xhat_hist = zeros(N,4);
    xhat_hist(1,:) = xhat';

    % ---- loop ----
    for k = 2:N
        % prediction
        xhat = F*xhat + B_u*a_imu(k-1);
        P    = F*P*F' + Q;

        % GPS position update
        if gps_new(k)
            z = p_gps(k);
            y  = z - H_gps*xhat;
            S  = H_gps*P*H_gps' + R_gps;
            K  = (P*H_gps')/S;
            xhat = xhat + K*y;
            P    = (I - K*H_gps)*P;
        end

        % GPS velocity update (filtered)
        if gps_new(k) && ~isnan(v_gps_f(k))
            z = v_gps_f(k);
            y  = z - H_gps_v*xhat;
            S  = H_gps_v*P*H_gps_v' + R_gps_v;
            K  = (P*H_gps_v')/S;
            xhat = xhat + K*y;
            P    = (I - K*H_gps_v)*P;
        end

        % Pitot update (EKF on q)
        if pitot_new(k)
            v_pred = xhat(2);
            w_pred = xhat(4);
            v_air_pred = v_pred - w_pred;

            if abs(v_air_pred) >= vmin_pitot
                z = q_pitot(k);

                h = 0.5*rho*(v_air_pred^2);
                H = [0, rho*v_air_pred, 0, -rho*v_air_pred];

                y = z - h;

                R_eff = R_pitot * (1 + (v0/max(abs(v_air_pred),0.2))^2);

                S = H*P*H' + R_eff;
                K = (P*H')/S;

                xhat = xhat + K*y;
                P    = (I - K*H)*P;
            end
        end

        xhat_hist(k,:) = xhat';
    end

    % ---- MAE metrics ----
    p_est = xhat_hist(:,1);
    v_est = xhat_hist(:,2);
    b_est = xhat_hist(:,3);
    w_est = xhat_hist(:,4);

    v_air_est  = v_est - w_est;
    v_air_true = v_true - w_true;

    MAE_p    = mean(abs(p_est - p_true));
    MAE_v    = mean(abs(v_est - v_true));
    MAE_w    = mean(abs(w_est - w_true));
    MAE_b    = mean(abs(b_est - b_true));
    MAE_vair = mean(abs(v_air_est - v_air_true));
end
