function [spl] = cubicSpline(x, f)

    close all;

    x = [-fliplr(x) x];
    f = [fliplr(f) f];
    
    x = x(:); f = f(:);
    N = length(x);
    h = zeros(N-1, 1);
    for i = 1:N-1
        h(i) = x(i+1) - x(i);
    end

    A = zeros(N-1, 1);
    B = zeros(N, 1);
    C = zeros(N-1, 1);
    D = zeros(N-1, 1);

    L = zeros(N, N); % Coefficient matrix for system of equation for B
    L(1, 1) = 1; L(N, N) = 1;
    F = zeros(N, 1); % RHS for system of equation for B
    F(1) = 0; F(N) = 0;

    % Generating the coefficient matrix and RHS
    for i = 2:N-1
        L(i, i-1) = h(i-1);
        L(i, i) = 2 * (h(i-1) + h(i));
        L(i, i+1) = h(i);
        F(i) = 3 * (f(i+1) - f(i)) / h(i) + 3 * (f(i-1) - f(i)) / h(i-1);
    end
    
    % Solving the equation
    B = linsolve(L, F);

    % Calculating A, C and D
    for i = 2:N
        A(i-1) = (B(i) - B(i-1)) / (3*h(i-1));
        C(i-1) = (f(i) - f(i-1)) / h(i-1) - (B(i) + 2*B(i-1)) * h(i-1) / 3;
    end
    D = f(1:N-1);

    % x-range on which cubic spline will be calculated
    temp = [];
    % The predicted values will be stored here
    spln = [];

    % Calculating the actual cubic polynomial's value at each point,
    % one polynomial at at time
    for i = 1:N-1
        temp(1000*i - 999 : 1000*i) = linspace(x(i), x(i+1), 1000)'; % 1000 evaluation points pre interval
        spln(1000*i - 999 : 1000*i) = A(i).*((temp(1000*i - 999 : 1000*i) - x(i)).^3) + B(i).*((temp(1000*i - 999 : 1000*i) - x(i)).^2) + C(i).*(temp(1000*i - 999 : 1000*i) - x(i)) + D(i);
    end

    % Generating the plots and comparing with MATLAB's spline function
    figure(1);
        subplot(1,2,1);
            plot(x, f, 'bo');
            hold on;
            plot(temp, spln, 'r');
            title("Cubic Spline: Airy Disc");
            xlabel("Distance from Center"); ylabel("Normalized Intensity");
        subplot(1,2,2);
            plot(x, f, 'bo');
            hold on;
            plot(temp, spline(x, f, temp), 'g');
            title("Cubic Spline: Airy Disc");
            xlabel("Distance from Center"); ylabel("Normalized Intensity");
    
end