consts.m_secco = 10; % kg
consts.m_propellente = 5; %kg
consts.isp = 300; % s
consts.g0 = 9.81;
consts.burn_time = 1 % s

consts.diametro = 0.10;
consts.CD = 0.4;           
consts.rho = 1.225;   
consts.area = pi * (consts.diametro/2)^2; 

mass = m_propellente + m_secco;
consts.m_dot = consts.m_propellente / consts.burn_time;

init_position = [0, 0, 0];
init_velocity = [0, 0, 0];

zenit = 1; % rad angolo tra l'asse del razzo, e l'asse Z
azimut = 0; % angolo tra l'asse Y, e la proiezione dell'asse del razzo sul piano XY
init_angle = [zenit, azimut]; 


%state vector:
% (1:3) posizione
% (4:6) velocita
% (7) massa
% (8, 9) angolo zenit, angolo azimut
state0 = [init_position' ; init_velocity'; mass; init_angle'];

% dm/m = v_e dv
% v_e m'/m = a 
options = odeset('Events', @(t, s) ground_event(t, s));
[t, state] = ode45(@(t, state) solver(t, state, consts), [0 100], state0, options)

figure(1);
plot(t, state(:, 7), 'LineWidth', 1.5);
grid on;
xlim([0 5])
xlabel('Tempo (s)');
ylabel('Massa (kg)');
title('Variazione della massa nel tempo');

figure(2)
plot3(state(:, 1)', state(:, 2)', state(:, 3)', '-o', 'LineWidth', 1.5, 'MarkerSize', 6);

xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

velocity_magnitude = sqrt(state(:, 4).^2 + state(:, 5).^2 + state(:, 6).^2);

figure(3);
plot(state(:, 7), velocity_magnitude, 'LineWidth', 1.5);
grid on;
xlabel('Massa (kg)');
ylabel('Velocità (m/s)');
title('Velocità in funzione della massa');
% state vettore colonna [x, y, x, vx, vy, vz]
% ZXY: Z: rivolto verso l'alto, X rivolto a destra, Y entrante 
function state_updated = solver(t, state, consts)
    vx = state(4);
    vy = state(5);
    vz = state(6);
    V = sqrt(vx^2 + vy^2 + vz^2);
    % acc - asse z
    gravity = -consts.g0;
    thrust_z = consts.isp * consts.g0 * consts.m_dot/state(7) * cos(state(8));
    
    % acc - asse x
    thrust_x = consts.isp * consts.g0 * consts.m_dot/state(7) * sin(state(8))*cos(state(9));

    thrust_y = consts.isp * consts.g0 * consts.m_dot/state(7) * sin(state(8))*sin(state(9));
    % Propellant is exhausted after the burn; thrust and mass flow then stop.
    if t >= consts.burn_time || state(7) <= consts.m_secco
        thrust_x = 0;
        thrust_y = 0;
        thrust_z = 0;
        consts.m_dot = 0;
    end
    if V > 1e-3
        D = 0.5 * consts.rho * V^2 * consts.CD * consts.area;
        % Le componenti dell'accelerazione di drag agiscono nella direzione opposta a [vx, vy, vz]
        drag_x = -(D / state(7)) * (vx / V);
        drag_y = -(D / state(7)) * (vy / V);
        drag_z = -(D / state(7)) * (vz / V);
    else
        drag_x = 0;
        drag_y = 0;
        drag_z = 0;
    end
    
    acc_x = thrust_x + drag_x;
    acc_y = thrust_y + drag_y;
    acc_z = thrust_z + gravity + drag_z;

    acceleration = [acc_x; acc_y; acc_z];
  
    state_updated = [state(4:6); acceleration; -consts.m_dot; 0; 0];
end

function [value, isterminal, direction] = ground_event(~, state)
value = state(3);       % Intercetta quando la coordinata Z (quota) è 0
isterminal = 1;         % 1 = ferma l'integrazione dell'ODE
direction = -1;         % -1 = rileva solo quando Z sta diminuendo (fase di discesa)
end