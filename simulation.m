clc
close all
clear all

%% Simulation Params

% Enviroment - UMa-NLoS - number 4 sorted
env = 'UMa_NLoS';
% Frequency (GHz)
freq_ghz = 3;
c = 299792458; % m/s
lambda = c / (freq_ghz * 1e9);
% Multipath components
N = 100;
% Elevation Angle
phi_bar = pi/4;
% RX Velocity
v = 5; % m/s


%% Delay Spread RMS
% Stats
[delay_mean, delay_std] = get_delay_stats(freq_ghz, env);
% Gaussian - only one sample to represent
delay_spread_rms_log = normrnd(delay_mean, delay_std, [1,1]);
delay_spread = 10^(delay_spread_rms_log);
delay_spread = 0.00000088316;
delay_spread_ns = 1e9*delay_spread

%% Generate Delay Samples
% Proportionality Factor
r_tau = get_prop_factor_delay(env);
mu_tau = delay_spread * r_tau;
% Delay samples
tau_n_absolute = exprnd(mu_tau, [N,1]);
tau_n_normalized = tau_n_absolute - min(tau_n_absolute);
tau_n = sort(tau_n_normalized);
% tau_n_ns = tau_n(1:5)*1e9

%% Generate Power Samples
% Power stats
[shadow_std, rice_mean, rice_std] = get_power_stats(env);
% Shadowing samples
shadow_samples = normrnd(0, shadow_std, [N,1]);
% Power samples
alpha_n_2 = exp(-tau_n .* (r_tau-1) ./ (r_tau*delay_spread)) .* 10.^(-shadow_samples ./ 10);
% alpha_n_2(1:5)

%% Rice factor
kr = 0;
kr_db = -inf;
if env == "UMi_LoS" || env == "UMa_LoS"
    kr_db = normrnd(rice_mean, rice_std, [1,1]);
    kr = 10 ^ (kr_db / 10);
end
kr

% Normalizing power
hat_omega_c = sum(alpha_n_2(2:N));
alpha_n_2 = (1/(kr + 1)) .* (alpha_n_2 ./ hat_omega_c);
alpha_n_2(1) = kr / (kr + 1);

% Verifying power = 1
omega_c = sum(alpha_n_2)
% Verifying Rice Factor
kr_check = alpha_n_2(1) ./ sum(alpha_n_2(2:N))
% Verifying Delay Spread
tau_bar = (1 / omega_c) * sum(tau_n .* alpha_n_2);
delay_spread
sigma_tau_check = sqrt((1 / omega_c) * sum(alpha_n_2 .* ((tau_n - tau_bar) .^ 2)))

% ------------------------------------------------
%% Plot Power Delay
% ------------------------------------------------
figure(1)
stem( tau_n / 1e-6, alpha_n_2, '^', 'Color', 'k', 'Linewidth', 1.5 );
hold on
stem( tau_n(1) / 1e-6, alpha_n_2(1), '^', 'Color', 'blue', 'Linewidth', 1.5 );

ylim([1e-10, 1]);
% xlim([0, 10]);

xlabel( 'Atraso multipercurso -- $\tau$ ($\mu$s)', 'Interpreter', 'Latex', 'Fontsize', 13 );
ylabel( 'Pot{\^{e}}ncia multipercurso', 'Interpreter', 'Latex', 'Fontsize', 13 );
grid on

ax = gca;
ax.TickLabelInterpreter = 'Latex';
ax.FontSize = 14;

% Put y-scale on log
set(ax,'yscal','log');


% --------------------------------------
%% Angles Stats
% --------------------------------------
[azimuth_mean, azimuth_std, elevation_mean, elevation_std] = get_angle_stats(freq_ghz, env);
% Azimuth Angles
sigma_theta_degree = 10 ^ (normrnd(azimuth_mean, azimuth_std, [1,1]));
sigma_theta_rad = sigma_theta_degree * (pi / 180);
max_alpha = max(alpha_n_2);
theta_n = 1.42 * sigma_theta_rad * sqrt(-log(alpha_n_2 ./ max_alpha));
% Adjusts
un = randsample([-1,1], N, true)';
% un(1:5)
yn = normrnd(0, sigma_theta_rad/7, [N, 1]);
% yn(1:5)
% Finally
theta_n = (un .* theta_n) + yn;
if env == "UMi_LoS" || env == "UMa_LoS"
    theta_n = theta_n - theta_n(1);
else
    theta_n(1) = 0;
end
% theta_n(1) = 0;
% theta_1_5 = theta_n(1:5)

% ------------------------------------------------
%% Plot Azimuth Angles
% ------------------------------------------------
figure(2)

if env == "UMi_NLoS" || env == "UMa_NLoS"
    min_alpha_n_2 = min( alpha_n_2(2:N) );
else
    min_alpha_n_2 = min( alpha_n_2 );
end

% min_alpha_n_2

for n = 1 : length( theta_n )
    power_x = 10 * log10( alpha_n_2(n) / min_alpha_n_2 );
    if n == 1 
        polarplot( [theta_n(n),theta_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'blue' );
    else
        polarplot( [theta_n(n),theta_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'black' );
    end
    hold on
end
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 2;
ax.ThetaTickLabel = {'0'; '30'; '60'; '90'; '120'; '150'; '180'; '210'; '240'; '270'; '300'; '330';};
title('Espalhamento angular em azimute', 'Interpreter', 'Latex');
grid on

% Plot Power with Azimuth Angles
figure(3)
stem( rad2deg( theta_n ), alpha_n_2, '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
hold on
stem( rad2deg( theta_n(1) ), alpha_n_2(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue' );
% ylim([1e-8,1]);
% xlim([-6 * (10 .^ sigma_theta_degree), 6 * (10 .^ sigma_theta_degree)]);
set(gca,'yscal','log');
xlabel('{\^{A}}ngulos de chegada em azimute ($^{\circ}$)', 'Interpreter', 'Latex');
ylabel('Pot{\^{e}}ncia', 'Interpreter', 'Latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
grid on

%% Elevation Angles
sigma_phi_degree = 10 ^ (normrnd(elevation_mean, elevation_std, [1,1]));
sigma_phi_rad = sigma_phi_degree * (pi / 180);
max_alpha = max(alpha_n_2);
phi_n = -sigma_phi_rad * log(alpha_n_2 ./ max_alpha);
% Adjusts
un = randsample([-1,1], N, true)';
% un(1:5)
yn = normrnd(0, sigma_phi_rad/7, [N, 1]);
% yn(1:5)
% Finally
phi_n = (un .* phi_n) + yn;
if env == "UMi_LoS" || env == "UMa_LoS"
    phi_n = phi_n - phi_n(1) + phi_bar;
else
    phi_n(1) = 0;
end
% phi_n(1) = 0;
% phi_1_5 = phi_n(1:5)

% ------------------------------------------------
%% Plot Elevation Angles
% ------------------------------------------------
figure(4)

if env == "UMi_NLoS" || env == "UMa_NLoS"
    min_alpha_n_2 = min( alpha_n_2(2:N) );
else
    min_alpha_n_2 = min( alpha_n_2 );
end

for n = 1 : length( phi_n )
    power_x = 10 * log10( alpha_n_2(n) / min_alpha_n_2 );
    if n == 1 
        polarplot( [phi_n(n),phi_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'blue' );
    else
        polarplot( [phi_n(n),phi_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'black' );
    end
    hold on
end
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 2;
ax.ThetaTickLabel = {'0'; '30'; '60'; '90'; '120'; '150'; '180'; '210'; '240'; '270'; '300'; '330';};
title('Espalhamento angular em elevacao', 'Interpreter', 'Latex');
grid on

% Plot Power with Azimuth Angles
figure(5)
stem( rad2deg( phi_n ), alpha_n_2, '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
hold on
stem( rad2deg( phi_n(1) ), alpha_n_2(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue' );
% ylim([1e-8,1]);
% xlim([-6 * (10 .^ sigma_theta_degree), 6 * (10 .^ sigma_theta_degree)]);
set(gca,'yscal','log');
xlabel('{\^{A}}ngulos de chegada em elevacao ($^{\circ}$)', 'Interpreter', 'Latex');
ylabel('Pot{\^{e}}ncia', 'Interpreter', 'Latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
grid on

%---------------------------------------
%% Directions
%---------------------------------------
rn = [cos(theta_n) .* sin(phi_n), sin(theta_n) .* sin(phi_n), cos(phi_n)];

if env == "UMi_NLoS" || env == "UMa_NLoS"
    rn(1, :) = 0;
end

% RX Direction Vector
theta_v = pi/4; % azimuth
phi_v = 0; % elevation

vrx = [cos(theta_v) .* sin(phi_v), sin(theta_v) .* sin(phi_v), cos(phi_v)];

%% Plot angle vectors
figure(6)
for i = 1 : N
    if i == 1
        style = "b-^";
    else
        style = "k-o";
    end
    plot3( [0 rn(i, 1)], [0 rn(i, 2)], [0 rn(i, 3)], style, 'LineWidth', 1.5 );
    hold on;
end
xlabel('x'), ylabel('y'), zlabel('z')
set(gca,'CameraPosition',[1 2 3]);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
grid on
style = "g->";
plot3( [0 vrx(1, 1)], [0 vrx(1, 2)], [0 vrx(1, 3)], style, 'LineWidth', 1.5);
blanck_labels = repmat({''}, 1, N-2);
legend(['LoS','NLoS', blanck_labels, 'V_{RX}'])

% -------------------------------------
%% Doppler Effects
% -------------------------------------
vn = (v / lambda) * (rn * vrx');
vn = sum(vn, 2);

% Plot doppler shift
figure(8)
stem( vn(1), alpha_n_2(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue');
hold on
stem( vn(2:N) , alpha_n_2(2:N), '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
% ylim([1e-8,1]);
% xlim([-6 * (10 .^ sigma_theta), 6 * (10 .^ sigma_theta)]);
set(gca,'yscal','log');
grid on
title("$v_{rx}="+num2str(v)+"$ m/s, $f_c = "+num2str(freq_ghz)+"$ GHz", 'Interpreter', 'Latex')
xlabel('Desvio Doppler - $\nu$ (Hz)', 'Interpreter', 'Latex');
ylabel('Potencia', 'Interpreter', 'Latex');
legend('LoS', 'NLoS')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
grid on

% --------------------------------------
%% Signal Analysis
% --------------------------------------
deltas_t = [1e-7, 1e-5, 1e-3];
n_samples = 1e5;


% Transmited Signal
t = zeros(n_samples, length(deltas_t));
signal_tx = zeros(n_samples, length(deltas_t));
for d = 1 : length(deltas_t)
    t(: , d) = linspace(0, 5*deltas_t(d), n_samples);
    [signal_tx(:, d)] = generate_pulse(0, t(: , d), deltas_t(d));
    
    % % Plot Transmited Signal
    % if d == 1 
    %     figure(9)
    % end
    % subplot(3, 1, d);
    % plot(t(:, d), abs(signal_tx(:, d)), 'Color', 'b', 'Linewidth', 1.5)
    % hold on
    % % xticks(-2*delta:delta:max(t))
    % % xticklabels({'$-2\delta_t$', '$-\delta_t$', '$0$', '$\delta_t$', '$2\delta_t$', '$3\delta_t$', '$4\delta_t$', '$5\delta_t$', '$6\delta_t$'})
    % title('Sinal Transmitido $\delta t = 10^{-7}s$', 'Interpreter', 'Latex')
    % ylabel('$s(t)$', 'Interpreter', 'Latex')
    % xlabel('$t$ [s]', 'Interpreter', 'Latex')
    % legend('TX')
    % ax = gca;
    % ax.TickLabelInterpreter = 'latex';
    % ax.FontSize = 14;
    % grid on

end

% --------------------------------------
%% Received Signal
% --------------------------------------


signal_rx = zeros(n_samples, N, length(deltas_t));
% signal_rx = zeros(N, length(deltas_t));
fig_init = 10;

for d = 1 : length(deltas_t)
    for i = 1 : N
        delayed_signal = generate_pulse(tau_n(i), t(:, d), deltas_t(d));
        phase_n = 2*pi * (((freq_ghz * 1e9 + vn(n)) * tau_n(n)) - (vn(n) * t(:, d)));
        signal_rx(:, i, d) = alpha_n_2(i) * exp(-phase_n * 1j) .* delayed_signal;
    end

    scattered_signal_rx(:, d) = sum(signal_rx(:, :, d), 2);
    abs_scat_signal_rx(:, d) = abs(scattered_signal_rx(:, d));
    normalized_signal_rx(:, d) = (abs_scat_signal_rx(:, d) - min(abs_scat_signal_rx(:, d))) / (max(abs_scat_signal_rx(:, d)) - min(abs_scat_signal_rx(:, d)));

    figure(fig_init)
    plot(t(:, d), abs(signal_tx(:, d)), 'Color', 'b', 'Linewidth', 1.5)
    hold on
    plot(t(:, d), normalized_signal_rx(:, d), 'Color', 'r', 'Linewidth', 1.2)
    xlim([0, 5*deltas_t(d)])
    % xticklabels({'$-2\delta_t$', '$-\delta_t$', '$0$', '$\delta_t$', '$2\delta_t$', '$3\delta_t$', '$4\delta_t$', '$5\delta_t$', '$6\delta_t$'})
    ylabel('$s(t)$', 'Interpreter', 'Latex')
    xlabel('$t$ [s]', 'Interpreter', 'Latex')
    title("$\delta t="+num2str(deltas_t(d)*1e6)+" \mu s$, $\sigma_{\tau}="+num2str(delay_spread_ns)+"\eta s$, Rice$="+num2str(kr_db)+"$dB", 'Interpreter', 'Latex')
    legend('TX', 'RX')
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 14;
    grid on

    fig_init = fig_init + 1;
end

%% Coherency bandwidth and coherency time
n_samples = 1e4;
k = logspace(-3, 10, n_samples);
t = logspace(-6, 0, n_samples);

omega_c = sum( alpha_n_2 );
autocorr_freq_n = (1 / omega_c) .* alpha_n_2 .* exp( -2i * pi * tau_n .* k );
autocorr_freq = sum(autocorr_freq_n, 1);

autocorr_time_n = (1 / omega_c) .* alpha_n_2 .* exp( 2i * pi * vn .* t);
autocorr_time = sum(autocorr_time_n, 1);

figure(13)
semilogx(k, abs(autocorr_freq), 'Color', 'r', 'Linewidth', 1.5)
hold on
yline(0.95, '-.k', 'Linewidth', 1.5)
yline(0.9, '--k', 'Linewidth', 1.5)
idx_095_freq = find(abs(autocorr_freq) < 0.95, 1,'first');
idx_090_freq = find(abs(autocorr_freq) < 0.9, 1, 'first');
xline(k(idx_095_freq), '-.k', 'Linewidth', 1.5)
xline(k(idx_090_freq), '--k', 'Linewidth', 1.5)
xlim([min(k), max(k)]);
ylabel('$\vert\rho_{TT}(\kappa, 0)\vert$', 'Interpreter', 'Latex')
xlabel('$\kappa$ [Hz]', 'Interpreter', 'Latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
title("$B_C(0,95)="+sprintf('%.2f', k(idx_095_freq)*1e-6)+"$MHz, $B_C(0,9)="+sprintf('%.2f', k(idx_090_freq)*1e-6)+"$MHz, $\sigma_{\tau}="+num2str(delay_spread_ns)+"\eta s$, Rice$="+num2str(kr_db)+"$dB", 'Interpreter', 'Latex', 'Fontsize', 10)
grid on

figure(14)
semilogx(t, abs(autocorr_time), 'Color', 'r', 'Linewidth', 1.5)
hold on
yline(0.95, '-.k', 'Linewidth', 1.5)
yline(0.9, '--k', 'Linewidth', 1.5)
idx_095_time = find(abs(autocorr_time) < 0.95, 1, 'first');
idx_090_time = find(abs(autocorr_time) < 0.9, 1, 'first');
xline(t(idx_095_time), '-.k', 'Linewidth', 1.5)
xline(t(idx_090_time), '--k', 'Linewidth', 1.5)
xlim([min(t), max(t)]);
ylabel('$\vert\rho_{TT}(0, \sigma)\vert$', 'Interpreter', 'Latex')
xlabel('$\sigma$ [s]', 'Interpreter', 'Latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
title("$T_C(0,95)="+sprintf('%.2f', t(idx_095_time)*1e3)+"$ms, $T_C(0,9)="+sprintf('%.2f', t(idx_090_time)*1e3)+"$ms, $\sigma_{\tau}="+num2str(delay_spread_ns)+"\eta s$, Rice$="+num2str(kr_db)+"$dB", 'Interpreter', 'Latex', 'Fontsize', 10)
grid on

% --------------------------------------
%% Functions
% --------------------------------------
function [delay_mean, delay_std] = get_delay_stats(freq_ghz, env)
    switch env
        case 'UMi_LoS'
            delay_mean = -0.24 * log10( 1 + freq_ghz ) - 7.14;
            delay_std = 0.38;

        case 'UMi_NLoS'
            delay_mean = -0.24 * log10( 1 + freq_ghz ) - 6.83;
            delay_std = -0.16 * log10( 1 + freq_ghz ) + 0.28;
            
        case 'UMa_LoS'
            delay_mean = -0.0963 * log10( 1 + freq_ghz ) - 6.955;
            delay_std = 0.66;
            
        case 'UMa_NLoS'
            delay_mean = -0.204 * log10( 1 + freq_ghz ) - 6.28;
            delay_std = 0.39;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [prop_factor_delay] = get_prop_factor_delay(env)
    switch env
        case 'UMi_LoS'
            prop_factor_delay = 3;

        case 'UMi_NLoS'
            prop_factor_delay = 2.1;
            
        case 'UMa_LoS'
            prop_factor_delay = 2.5;
            
        case 'UMa_NLoS'
            prop_factor_delay = 2.3;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [shadow_std, rice_mean, rice_std] = get_power_stats(env)
    switch env
        case 'UMi_LoS'
            shadow_std = 4;
            rice_mean = 9;
            rice_std = 5;
            
        case 'UMi_NLoS'
            shadow_std = 7.82;
            rice_mean = -inf;
            rice_std = -inf;
            
        case 'UMa_LoS'
            shadow_std = 4;
            rice_mean = 9;
            rice_std = 3.5;
            
        case 'UMa_NLoS'
            shadow_std = 6;
            rice_mean = -inf;
            rice_std = -inf;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [azimuth_mean, azimuth_std, elevation_mean, elevation_std] = get_angle_stats(freq_ghz, env)
    switch env
        case 'UMi_LoS'
            azimuth_mean = -0.08 * log10( 1 + freq_ghz ) + 1.73;
            azimuth_std = 0.014 * log10( 1 + freq_ghz ) + 0.28;
            elevation_mean = -0.1 * log10( 1 + freq_ghz ) + 0.73;
            elevation_std = -0.04 * log10( 1 + freq_ghz ) + 0.34;

        case 'UMi_NLoS'
            azimuth_mean = -0.08 * log10( 1 + freq_ghz ) + 1.81;
            azimuth_std = 0.05 * log10( 1 + freq_ghz ) + 0.3;
            elevation_mean = -0.04 * log10( 1 + freq_ghz ) + 0.92;
            elevation_std = -0.07 * log10( 1 + freq_ghz ) + 0.41;
            
        case 'UMa_LoS'
            azimuth_mean = 1.81;
            azimuth_std = 0.2;
            elevation_mean = 0.95;
            elevation_std = 0.16;
            
        case 'UMa_NLoS'
            azimuth_mean = -0.27 * log10( freq_ghz ) + 2.08;
            azimuth_std = 0.11;
            elevation_mean = -0.3236 * log10( freq_ghz ) + 1.512;
            elevation_std = 0.16;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [signal] = generate_pulse( delay, t, pulse_width )
    signal = zeros(length(t), 1);
    for i = 1 : length(t)
        if t(i) >= delay && t(i) <= pulse_width + delay
            signal(i) = 1;
        end
    end
end