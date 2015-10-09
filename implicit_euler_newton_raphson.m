% =========================================================================
% Solve M x_ddot + B x_dot + K x = F
% Implicit Euler with Newton Raphson
% =========================================================================

function implicit_euler_newton_raphson()
    clear all;
    close all;
    clc;
    
    figure;
    hold;
    
    K = 100;           % spring stiffness. N/m
    M = 5;             % mass, kg
    B = 2*sqrt(K*M);  %critical damping
    %B = 0;   
    
    % initial condition:
    x0 = 0.1;
    v0 = 0;
    F0 = 0;
    
    h = 0.2;
    t_start = 0;
    t_end = 10;
    t=t_start:h:t_end;
    
    % Use ode45, 1kHz as ground truth:    
    [tode45,xode45]=ode45(@msd, [t_start:0.001:t_end], [x0 v0], [], ...
                          M, B, K);
    plot(tode45,xode45(:,1), '--r')
    
    for k = 1 : length(t)
        % Newton-Raphosn relies on good intial value
        % As an initial guess, a 1-step forward Euler is used
        %v1 = v1_hat(M, B, K, 0, x0, v0, h);     

        % Initializing v1 = 0 also works. It just needs more iterations
        v1 = 0;        
        
        % Newton-Rapshon
        % x(k+1) = x(k) - g(x(k))/g'(x(k))
        % Since the system is linear, g'(x(k)) = constant
        g_d =  M/h + B + K * h;
        while(1) % Newton-Raphson iteration
            v1_ = v1 - (g(M, B, K, F0, x0, v0, v1, h) / g_d);            
            if abs(v1 - v1_) <0.0001
                break;
            end
            v1 = v1_;
        end                        
        
        x1 = x0 + v0 * h;   
        x0 = x1;
        v0 = v1;        
        
        x(k, :) = [x1 v1];
    end
    plot(t, x(:,1), 'b')   
    
    legend('ode45', 'Implicit Euler');
    xlabel('Time (s)')
end

function output = v1_hat(M, B, K, F, x0, v0, h)
    v0_dot = (F - B * v0 - K * x0) / M;
    output = v0 + h * v0_dot;
end

function output = g(M, B, K, F, x0, v0, v1, h)
    output = M*(v1-v0)/h-F+B*v1+K*(x0+h*v0);
end

function xdot=msd(t,x, M, B, K)
    xdot_1 = x(2);
    xdot_2 = -(B/M)*x(2) - (K/M)*x(1);

    xdot = [xdot_1 ; xdot_2 ];
end

