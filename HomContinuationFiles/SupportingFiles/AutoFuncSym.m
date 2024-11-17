function f = AutoFuncSym(x_in, p_in)
    % Generates a COCO-compatible symbolic function handle for a dynamical system.
    % Defines symbolic equations and parameters, then converts them into a COCO-compatible function.
    % Outputs:
    %   - f: A COCO-compatible function handle.

    %%% Symbolic variables
    syms y_N y_E x_B x_N          % State variables
    syms m eta k2 ep              % System parameters

    %%% System constants and parameters
    SV = 10^6;                    % Svedrups
    sig = 2.1 * 10^4 * SV;        % Advection strength
    alphaT = 1.7 * 10^(-4);       % Thermal expansion coefficient
    alphaS = 0.8 * 10^(-3);       % Salinity contraction coefficient
    TaN = 7;                      % Reference temperature (north)
    TaL = 25;                     % Reference temperature (low)
    T0 = 2.65;                    % Reference base temperature
    S0 = 35;                      % Reference salinity
    SB = 34.538;                  % Salinity in region B
    y_B = alphaS / alphaT * (SB - S0) / (TaN - T0);  % Computed salinity-related parameter
    VN = 7.2106 * 10^15;          % Volume for region N
    VL = 6.3515 * 10^16;          % Volume for region L
    VB = 18 * VN;                 % Volume for region B
    g = 3.17 * 10^-8;             % Gravitational constant (scaled)
    phi = alphaT * (TaN - T0);    % Temperature-dependent factor
    x_E = (TaL - T0) / (TaN - T0);% Temperature-dependent variable
    TildeVE = VL / VN;            % Volume ratio between regions L and N
    TildeVB = VB / VN;            % Volume ratio between regions B and N
    W = 5.456 * SV / sig;         % Wind forcing term
    delta_N = g * VN / sig;       % Temperature exchange coeff.

    %%% Advective functions
    ps = phi * ((y_N - x_N) - (y_E - x_E)); % Advection term (temperature gradient-based)
    Psi = ps * tanh(ps / ep);              % Smoothed nonlinear advection

    %%% Convective functions
    k2n = k2;                              % Convective coefficient (region N)
    k2l = k2 / 2.0;                        % Convective coefficient (region L)
    k1 = 10^-2 * SV / sig;                 % Minimum diffusion coefficient
    % Convective terms with smoothing
    kappaN = k1 + 0.5 * (k2n - k1) * (1 + tanh(((y_N - y_B) - (x_N - x_B) - eta) / ep)); 
    kappaE = k1 + 0.5 * (k2l - k1) * (1 + tanh(((y_E - y_B) - (x_E - x_B) - eta) / ep)); 

    %%% System equations
    % y_N (North region temperature equation)
    ynp = m + ps / 2 * (y_E - y_B) + Psi / 2 * (y_E + y_B - 2 * y_N) + ...
          W * (y_E - y_N) - kappaN * (y_N - y_B);

    % y_E (Equatorial region temperature equation)
    ylp = 1 / TildeVE * (-m + ps / 2 * (y_B - y_N) + Psi / 2 * (y_B + y_N - 2 * y_E) - ...
          W * (y_E - y_N) - kappaE * (y_E - y_B));

    % x_B (Salinity in region B)
    xbp = 1 / TildeVB * (ps / 2 * (x_N - x_E) + Psi / 2 * (x_N + x_E - 2 * x_B) + ...
          kappaN * (x_N - x_B) + kappaE * (x_E - x_B));

    % x_N (Salinity in region N)
    xnp = -delta_N * (x_N - 1) + ps / 2 * (x_E - x_B) + Psi / 2 * (x_E + x_B - 2 * x_N) + ...
          W * (x_E - x_N) - kappaN * (x_N - x_B);

    %%% COCO-compatible output
    % Convert symbolic equations to numerical functions with COCO compatibility
    f = sco_sym2funcs([ynp; ylp; xbp; xnp]*10^2, ...  % System equations scaled
                      {[y_N; y_E; x_B; x_N], [m; eta; k2; ep]}, ... % State and parameters
                      {'x', 'p'}, ...                  % Labels for states and parameters
                      'filename', 'sys', ...           % Save to 'sys.m'
                      'maxorder', 3);                  % Derivatives up to 3rd order
end
