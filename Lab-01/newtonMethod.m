%===================== Newton-Raphson Method Root Finder =====================%
%  
%  AUTHOR: Prakhar Srivastava -2016B5A70438G- BITS Pilani K.K. Birla Goa Campus
%
%  USAGE: [root, iter, tol] = newtonMethod(initialGuess, tolerance, func, funcPrime)
%
%  Given a function "func" and its derivative "funcPrime", approximately 
%  determines the root of "func" with relative tolerance "tolerance".
%  Also takes the parameter "initialGuess" to start the root-searching from,
%  and return actual tolerance as "tol".
%
%=============================================================================%

function [root, iter, tol] = newtonMethod(initialGuess, tolerance, func, funcPrime)

    x0 = initialGuess;
    maxIter = 10e6;                % Prevent infinite loop
    iter = 0;
    tol = 0;

    while iter < maxIter

        y0 = func(x0);              % Calculate function's value
        yPrime = funcPrime(x0);     % Calculate derivative's value

        x1 = x0 - y0 ./ yPrime;      % Move closer to root
        tol = abs(x0 - x1) ./ x0;    % Current relative tolerance
        iter += 1;

        if tol < tolerance % If desired tolerance has been achieved
            break;
        endif

        x0 = x1;                    % Prepare for next loop

    endwhile

    root = x1;                      % Return this final value

endfunction