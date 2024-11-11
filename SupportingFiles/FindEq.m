function [roots] = FindEq(mu, eta, numRoots)
    % Computes roots of the UNFORCED system
    % Args:
    %   - mu and eta : parameters.
    %   - numRoots (optional) : number of roots returned. Defaults to one. 
    %

    % Default system parameters 
    SV=10^6;              
    sig=2.1*10^4*SV; 

    % Number of guesses + root array
    numGuesses = 20;
    roots = [];

    % Root guesses /w scaling
    L = 15.0;
    x1_range = L*(rand(1, numGuesses)); 
    x2_range = L*rand(1, numGuesses); 
    x3_range = L*rand(1, numGuesses);

    % Loop over the guesses to try converge to a root
    for i = 1:numGuesses
        % options = optimoptions('fsolve', 'Display','off', 'Algorithm', 'levenberg-marquardt', 'InitDamping', 0.001, 'FunctionTolerance', 10^-10);                
        % [sol, ~, exitflag, ~] = fsolve(@(x) AutoFuncODE(x, [mu, eta, 8.0*SV/sig]), [x1_range(i), x2_range(i), x3_range(i)], options);

        % Options for the Trust-Region Reflective algorithm
        options = optimoptions('fsolve', ...
            'Algorithm', 'trust-region-reflective', ...
            'Display', 'off', ...
            'FunctionTolerance', 1e-14, ...
            'StepTolerance', 1e-16, ...
            'OptimalityTolerance', 1e-14, ...
            'MaxFunctionEvaluations', 1e5, ...
            'MaxIterations', 1e2);
        
        % Solve using fsolve with enhanced accuracy settings
        [sol, ~, exitflag, ~] = fsolve(@(x) AutoFuncODE(x, [mu, eta, 8.0*SV/sig, 0.015]), ...
            [x1_range(i), x2_range(i), x3_range(i), x3_range(i)], options);

        % fsolve converged
        if exitflag > 0 % 
            roots = [roots; sol(1) sol(2) sol(3) sol(4)];
        end
    end
    
    % Only return the unique roots to `roundThresh' decimal places
    roundThresh = 5;
    roots = unique(round(roots, roundThresh), 'rows')

    % Default the number of roots returned to ALL. 
    if nargin < 3
        numRoots = length(roots);
    end

    % Catch the case when `numRoots' is too large for `roots'
    try
        roots = roots(1:end,:);
    catch
        fprintf("Number of roots chosen invalid. Returned one root instead.")
        roots = roots(1,:);
    end
end