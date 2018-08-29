% ============= Driver script for Newton-Raphson Method Root Finder =============
%
%  AUTHOR: Prakhar Srivastava -2016B5A70438G- BITS Pilani K.K. Birla Goa Campus
%
%  GUIDE: The script has been setup to solve the transcendental equation 
%  obtained from the Kronig Penney model. It will print the first five roots for
%  k = 0 and k = pi/2. Then it will generate a plot of the values of the first
%  five energy eigenstates vs. k.
%  User may tweak the ranges of parameters k and alpha.
%
%  NOTE: Script running time depends linearly on size of k. So do not set k
%  very large in size.
%
%================================================================================

tic; clc; clear all; close all;

% k can be set from here
k = [0 : 0.075 : 10*pi];
% location of pi/2 in k. Needs to be changed if k is changed.
loc = 1 + round(pi/(2 * 0.075));

% Various starting points;
alpha = [0.01 : 0.01 : 30];

% Some necessary variables
root = zeros(length(k), 5);
tolerance = 10e-8;
m_e = 9.e-31;
hbar = 1.0545718e-34;

for i = 1:length(k)

    % Target function and it's derivative
    f = @(x) (cos(k(i)) - cos(x) - 50 * sin(x) ./ x);
    fPrime = @(x) (sin(x) - 50 * (cos(x) ./ x - sin(x) ./ (x.^2)));

    % Evaluating roots by starting at various points
    % Can add additional return parameters if needed (check usage in documentation)
    [Roots] = newtonMethod(alpha, tolerance, f, fPrime);

    % Uncomment to display all roots for various starting points
    % fprintf("\n\nFor k = %f\n\n\[Initial Guess\]\t\[Root Found\]\n", k(i));
    % disp([alpha' Roots']);

    % Extract first five roots for various starting points
    root(i, :) = (unique(abs(Roots)))(1:5);
    
    % Uncomment to see the first five roots
    % fprintf("\n\nFor k = %f, first five roots are:\n", k(i));
    % disp(root(i, :));

endfor

% Print out the first five roots for k = 0 & pi/2
fprintf("For k = 0, first five roots are:\n");
disp(root(1, :));
fprintf("\nFor k = pi/2, first five roots are:\n");
disp(root(loc, :));

% Plot Energy vs. k
fprintf("\nGenerating the plot for Energy vs. k...\n\n");
plot(k, (root .* hbar).^2 / (2 * m_e), '.');
xlabel("K"); ylabel("First 5 allowed states");

toc;