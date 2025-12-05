% Wing Position Optimization
% =========================================


clc;
clear;
close all;

%% Optimization Setup
x0 = [1850, 720];
lb = [1600, 650];
ub = [2100, 775];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'UseParallel', false, 'SpecifyObjectiveGradient', false);

%% Run Optimization
[x_opt, T_opt, exitflag, output] = fmincon(@lap_time_wrapper, x0, [], [], [], [], lb, ub, @wing_constraints, options);

x_s = x_opt(1);
z_s = x_opt(2);

%% Display Results
fprintf('Exit flag: %d\n', exitflag);
fprintf('Optimal wing fore-aft x (mm): %.4f\n', x_s);
fprintf('Optimal wing height z   (mm): %.4f\n', z_s);
fprintf('Optimal lap time T(x,z) (s): %.6f\n', T_opt);
fprintf('========================================\n');

%% Contour Plot
nx = 80; nz = 80;
x_vec = linspace(lb(1), ub(1), nx);
z_vec = linspace(lb(2), ub(2), nz);
[X, Z] = meshgrid(x_vec, z_vec);

T_grid  = lap_time(X, Z);
LC_grid = load_cap(X, Z);

figure('Position', [100, 100, 800, 600]);
contour(X, Z, T_grid, 30, 'ShowText', 'on');
hold on;
contour(X, Z, LC_grid - 0.70, [0 0], 'LineWidth', 2, 'LineColor', 'k');
plot(x_s, z_s, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Wing fore-aft x (mm)');
ylabel('Wing height z (mm)');
title('Lap Time Contours with LC(x,z) = 0.70 and Optimal Point');
colorbar;
legend('Lap time', 'LC = 0.70', 'Optimal', 'Location', 'best');
grid on;
hold off;

%% Sensitivity Analysis - X Direction
npts = 30;
x_line = linspace(lb(1), ub(1), npts);
z_line_const = z_s * ones(1, npts);
T_x_line = lap_time(x_line, z_line_const);

p_x = polyfit(x_line, T_x_line, 2);
T_x_fit = polyval(p_x, x_line);

figure('Position', [120, 120, 800, 600]);
plot(x_line, T_x_line, 'o', 'MarkerSize', 6, 'DisplayName', 'Model samples');
hold on;
plot(x_line, T_x_fit, '-', 'LineWidth', 2, 'DisplayName', 'Quadratic fit');
xlabel('Wing fore-aft x (mm)');
ylabel('Lap time T(x,z_{opt}) (s)');
title(sprintf('Lap Time vs x at z = %.1f mm', z_s));
legend('Location', 'best');
grid on;
hold off;

%% Sensitivity Analysis - Z Direction
z_line = linspace(lb(2), ub(2), npts);
x_line_const = x_s * ones(1, npts);
T_z_line = lap_time(x_line_const, z_line);

p_z = polyfit(z_line, T_z_line, 2);
T_z_fit = polyval(p_z, z_line);

figure('Position', [140, 140, 800, 600]);
plot(z_line, T_z_line, 'o', 'MarkerSize', 6, 'DisplayName', 'Model samples');
hold on;
plot(z_line, T_z_fit, '-', 'LineWidth', 2, 'DisplayName', 'Quadratic fit');
xlabel('Wing height z (mm)');
ylabel('Lap time T(x_{opt},z) (s)');
title(sprintf('Lap Time vs z at x = %.1f mm', x_s));
legend('Location', 'best');
grid on;
hold off;

%% Helper Functions

function T_val = lap_time_wrapper(xz)
    % Wrapper for fmincon (expects vector input)
    T_val = lap_time(xz(1), xz(2));
end

function CL_val = CL_fun(x, z)
    % Lift coefficient as function of position
    CL_val = -(0.02 + 0.05 * exp(...
        -((x - 1850) / 160).^2 - ((z - 740) / 55).^2));
end

function CD_val = CD_fun(x, z)
    % Drag coefficient as function of position
    CD_val = 0.90 ...
        - 0.15 * exp(-((x - 1850) / 180).^2) ...
        - 0.08 * exp(-((z - 710) / 60).^2) ...
        + 1.5e-7 * (x - 2050).^2;
end

function T_val = lap_time(x, z)
    % Lap time model
    CL = CL_fun(x, z);
    CD = CD_fun(x, z);
    
    u = -CL;
    u_plus = max(u, 0);
    
    F_S = atan(8 * u_plus);
    F_D = CD;
    F_I = CD .* u;
    
    T_val = 60 - 18 * F_S + 7 * F_D + 0.6 * F_I;
end

function LC_val = load_cap(x, z)
    % Load capacity constraint
    CL = CL_fun(x, z);
    CD = CD_fun(x, z);
    
    u = -CL;
    u_plus = max(u, 0);
    
    LC_val = CD + 0.35 * u_plus;
end

function [c, ceq] = wing_constraints(xz)
    % Nonlinear constraints for optimization
    LC = load_cap(xz(1), xz(2));
    c = LC - 0.70;  % LC must be <= 0.70
    ceq = [];       % No equality constraints
end