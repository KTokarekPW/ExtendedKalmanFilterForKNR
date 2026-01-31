%% EKF: IMU (propagacja) + GPS (pozycja) + PITOT (cisnienie q, nieliniowy pomiar)
clc; clear; close all;

%% ========== 1) USTAWIENIA SYMULACJI ==========
dt_imu = 0.01;   % [s]
T      = 60;     % [s]
t      = (0:dt_imu:T)';
N      = numel(t);
acc    = 0.6;
Nmc = 10000;                 % liczba uruchomień (Monte Carlo)

MAE_a_all = zeros(Nmc,1);
MAE_v_all = zeros(Nmc,1);
MAE_p_all = zeros(Nmc,1);  % Wektory pod MAE do zrobienia średniej z 

dt_gps   = 0.2;
gps_step = round(dt_gps/dt_imu);

dt_pitot   = 0.005;
pitot_step = round(dt_pitot/dt_imu);

rho = 1.225;

%% ======= "PRAWDA" =======
a_true = zeros(N,1);
a_true(t>=2  & t<8)  =  1.2*acc;
a_true(t>=15 & t<18) = -1.5*acc;
a_true(t>=25 & t<35) =  0.6*acc;
a_true(t>=45 & t<50) = -0.8*acc;

v_true = zeros(N,1);
p_true = zeros(N,1);
for k = 2:N
    v_true(k) = v_true(k-1) + a_true(k-1)*dt_imu;
    p_true(k) = p_true(k-1) + v_true(k-1)*dt_imu + 0.5*a_true(k-1)*dt_imu^2;
end
parpool;              % start puli workerów
parfor mc = 1:Nmc
%% ======= 3) MODELE SZUMÓW =======
% IMU
sigma_a_meas = 0.20;   % [m/s^2]
sigma_b_rw   = 0.03;   % [m/s^2/sqrt(s)]

b_true = zeros(N,1);
for k = 2:N
    b_true(k) = b_true(k-1) + sigma_b_rw*sqrt(dt_imu)*randn;
end
a_imu = a_true + b_true + sigma_a_meas*randn(N,1);

% GPS
sigma_gps_pos = 1.0;   % [m]
p_gps = p_true + sigma_gps_pos*randn(N,1);

gps_new = false(N,1);
gps_new(1:gps_step:end) = true;

% PITOT: mierzymy q = 0.5*rho*v^2 + noise
dq_3sigma = 1.0;           % [Pa] (Twoje założenie "3-sigma")
sigma_q   = dq_3sigma/3;   % [Pa] 1-sigma do R

q_true  = 0.5*rho*v_true.^2;
q_pitot = q_true + sigma_q*randn(N,1);

pitot_new = false(N,1);    %Tworzy wektor boleanów, gdzie każda wartość wektora odnosi się do 
pitot_new(1:pitot_step:end) = true;

%% ======= 4) FILTR (stan: [p; v; b]) =======
F = [1 dt_imu -0.5*dt_imu^2;
     0 1     -dt_imu;
     0 0      1];

B_u = [0.5*dt_imu^2; dt_imu; 0];

H_gps = [1 0 0];      % GPS -> p

R_gps   = sigma_gps_pos^2;
R_pitot = sigma_q^2;  % w Pa^2 (stałe)

% Q
q_a = (0.50)^2;
q_b = (sigma_b_rw^2)*dt_imu;

Q = [0.25*dt_imu^4*q_a, 0.5*dt_imu^3*q_a, 0;
     0.5*dt_imu^3*q_a,  dt_imu^2*q_a,     0;
     0,                0,                q_b];

% Init
xhat = [p_gps(1); 0; 0];
P    = diag([25, 25, 0.01]);

xhat_hist = zeros(N,3);
xhat_hist(1,:) = xhat';

%% ======= 5) PETLA FILTRU =======
I = eye(3);

vmin_pitot = 1.5;  % opcjonalny próg, żeby nie bawić się w pitot przy v≈0

for k = 2:N

    % --- PREDYKCJA (IMU) ---
    xhat = F*xhat + B_u*a_imu(k-1);
    P    = F*P*F' + Q;

    % --- UPDATE GPS ---
    if gps_new(k)
        z = p_gps(k);

        y  = z - H_gps*xhat;
        S  = H_gps*P*H_gps' + R_gps;
        Kgain = (P*H_gps')/S;

        xhat = xhat + Kgain*y;
        P    = (I - Kgain*H_gps)*P;
    end

    % --- UPDATE PITOT (EKF) ---
    if pitot_new(k)
        v_pred = xhat(2);

        % opcjonalnie: ignoruj pitot przy bardzo małej predykcji v
        if abs(v_pred) >= vmin_pitot

            z = q_pitot(k);

            % h(x) = 0.5*rho*v^2
            h = 0.5*rho*v_pred^2;

            % Jacobian: dh/dx = [0, rho*v, 0]
            H = [0, rho*v_pred, 0];

            y  = z - h;
            S  = H*P*H' + R_pitot;
            Kgain = (P*H')/S;

            xhat = xhat + Kgain*y;
            P    = (I - Kgain*H)*P;
        end
    end

    xhat_hist(k,:) = xhat';
end

%% ======= 7) MAE =======
MAE_a = mean(abs(a_imu - a_true));
MAE_v = mean(abs(xhat_hist(:,2) - v_true));
MAE_p = mean(abs(xhat_hist(:,1) - p_true));

MAE_a_all(mc) = MAE_a;
MAE_v_all(mc) = MAE_v;
MAE_p_all(mc) = MAE_p;
end

fprintf('\n===== STATYSTYKI Z %d URUCHOMIEN =====\n', Nmc);

% v
fprintf('v  median(MAE): %.4f\n', median(MAE_v_all));
fprintf('v  IQR (25-75): %.4f .. %.4f\n', prctile(MAE_v_all,25), prctile(MAE_v_all,75));
fprintf('v  p95: %.4f\n', prctile(MAE_v_all,95));
fprintf('v  std: %.4f\n\n', std(MAE_v_all));

% p
fprintf('p  median(MAE): %.4f\n', median(MAE_p_all));
fprintf('p  IQR (25-75): %.4f .. %.4f\n', prctile(MAE_p_all,25), prctile(MAE_p_all,75));
fprintf('p  p95: %.4f\n', prctile(MAE_p_all,95));
fprintf('p  std: %.4f\n', std(MAE_p_all));