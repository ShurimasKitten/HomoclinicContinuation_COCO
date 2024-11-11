function [probSettings, thmEq, thmPO, thmHB, thmSN, thmHom, thmSNPst, thmSNPun, thmPDst] = loadDefaultSettings()
    % loadDefaultSettings initializes default settings for continuation analysis.
    %
    % Outputs:
    %   probSettings - Structure containing problem and continuation settings
    %   thmEq        - Plot theme for equilibria
    %   thmPO        - Plot theme for periodic orbits
    %   thmHB        - Plot theme for Hopf bifurcations
    %   thmSN        - Plot theme for saddle-node bifurcations
    %   thmHom       - Plot theme for homoclinic orbits
    %   thmSNPst     - Plot theme for stable saddle-node bifurcations
    %   thmSNPun     - Plot theme for unstable saddle-node bifurcations
    %   thmPDst      - Plot theme for stable period-doubling bifurcations

    % Initialize symbolic functions
    F = AutoFuncSym('');
    f = struct('f', F(''), ...
               'dfdxhan', F('x'), ...
               'dfdx',  F('x'), ...
               'dfdphan', F('p'),...
               'Dfdxdx', F({'x*v','x*v'}), ...
               'Dfdpdx', F({'p*v','x*v'}), ...
               'dfdx_dx', F('x*v'), ...
               'Dfdpdp', F({'p*v','p*v'}),...
               'dfdxdx_dx', F({'x','x*v'}),...
               'dfdpdx_dx', F({'p','x*v'}), ...
               'Dfdxdxhan', F({'x*v','x*v'}), ...
               'Dfdxdxdxhan', F({'x*v','x*v','x*v'}),...
               'nanflag', 1);

     
    % Initialize problem settings structure
    probSettings = struct;
    probSettings.labs = {};          % Labels
    probSettings.UZ = {};            % User-defined zones or similar
    probSettings.f = f;              % Function definitions
    probSettings.funcODE = @AutoFuncODE;    % Store ODE 

    %% Continuation Settings
    probSettings.contSettings = struct;
    probSettings.contSettings.norm = Inf;               % Normalization parameter
    probSettings.contSettings.PtMX = [500 500];         % Maximum points
    probSettings.contSettings.h0 = 1e-3;               % Initial step size
    probSettings.contSettings.h_min = 1e-8;            % Minimum step size
    probSettings.contSettings.h_max = 2e-3;            % Maximum step size
    probSettings.contSettings.NAdapt = 7;              % Adaptation steps
    probSettings.contSettings.contVar = {{'mu', 'eta'}, {[-1 1], [-15 5]}}; % Continuation variables and their ranges

    %% Correction Settings
    probSettings.corrSettings = struct;
    probSettings.corrSettings.TOL = 1e-6;              % Tolerance for correction
    probSettings.corrSettings.ItMX = 15;               % Maximum iterations

    %% Collocation Settings
    probSettings.collSettings.NTST = 50;                % Number of test points
    probSettings.collSettings.NCOL = 4;                 % Number of collocation points

    %% Periodic Orbit (PO) Settings
    probSettings.contSettings.PO = struct;
    probSettings.contSettings.PO.contVar = {{'mu', 'eta', 'po.period'},  [-1 1]}; % Continuation variables for PO

    %% Plot Settings
    % Equilibrium Plot Theme
    thmEq = ep_plot_theme('ep');
    thmEq.special = {'HB', 'SN', 'BP'}; % Special markers for Hopf, Saddle-Node, and Bifurcation Points
    thmEq.HB = {'kp', 'MarkerFaceColor', 'red', 'MarkerSize', 10};
    thmEq.SN = {'kp', 'MarkerFaceColor', 'black', 'MarkerSize', 10};
    thmEq.BP = {'diamond', 'MarkerFaceColor', 'blue', 'MarkerSize', 8};

    % Periodic Orbit Plot Theme
    thmPO = po_plot_theme('po');
    thmPO.special = {'SN', 'PD', 'BP', 'TR'}; % Special markers
    thmPO.TR = {'kp', 'MarkerFaceColor', 'cyan', 'MarkerSize', 10};
    thmPO.SN = {'kp', 'MarkerFaceColor', 'black', 'MarkerSize', 10};
    thmPO.PD = {'kp', 'MarkerFaceColor', [0.7176, 0.2745, 1.0000], 'MarkerSize', 10};
    thmPO.BP = {'diamond', 'MarkerFaceColor', [0.8000, 0.8000, 0.8000], 'MarkerSize', 10};
    thmPO.lspec = {{'b-', 'LineWidth', 1.5}, {'b-.', 'LineWidth', 1.5}};

    % Saddle-Node Bifurcation Plot Theme
    thmSN = ep_plot_theme('ep.SN');
    thmSN.special = {'HP', 'BP'};  % Special markers
    thmSN.HP = {'kp', 'MarkerFaceColor', 'red', 'MarkerSize', 8};
    thmSN.BP = {'kp', 'MarkerFaceColor', 'red', 'MarkerSize', 8};
    thmSN.lspec = {{'b-', 'LineWidth', 1.5, 'color', 'black'}, {'b-.', 'LineWidth', 1.5, 'color','black'}};

    % Stable Saddle-Node Bifurcation Plot Theme
    thmSNPst = struct;
    thmSNPst.special = {};
    thmSNPst.lspec = {'color', [0.6431, 0.2588, 0.7020], 'linestyle', '-', 'LineWidth', 1.5};

    % Unstable Saddle-Node Bifurcation Plot Theme
    thmSNPun = struct;
    thmSNPun.special = {};
    thmSNPun.lspec = {'color', [0.6745, 0.5059, 0.7020], 'linestyle', ':', 'LineWidth', 1.5};

    % Stable Period-Doubling Bifurcation Plot Theme
    thmPDst = struct;
    thmPDst.special = {};
    thmPDst.lspec = {'color', [0.3020, 0.7451, 0.9333], 'linestyle', '-', 'LineWidth', 1.5};

    % Hopf Bifurcation Plot Theme
    thmHB = ep_plot_theme('ep.HB');
    thmHB.special = {'GH', 'BT', 'UZ', 'HH'}; % Special markers
    thmHB.ustab = 'L1';
    thmHB.ustabfun = @(x) 1 + (~isnan(x) & x > 0) + 2 * isnan(x); % Update strategy
    thmHB.usept = {'GH','BT', 'BTP', 'BP'}; % Points to use
    thmHB.GH = {'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 8};
    thmHB.BT = {'square', 'MarkerFaceColor', 'k', 'MarkerSize', 8};
    thmHB.BTP = {'square', 'MarkerFaceColor', 'k', 'MarkerSize', 1};
    thmHB.BP = {'square', 'MarkerFaceColor', 'k', 'MarkerSize', 1};
    thmHB.HH = {'diamond', 'MarkerFaceColor', 'cyan', 'MarkerSize', 8};
    thmHB.UZ = {'diamond', 'MarkerFaceColor', 'black', 'MarkerSize', 10};
    thmHB.lspec = {{'r-', 'LineWidth', 1.5}, {'r-.', 'LineWidth', 1.5}};

    % Homoclinic Orbit Plot Theme
    red    = [0.7020, 0.0980, 0.0980];
    blue   = [0, 0.4471, 0.7412];
    purple = [0.7176, 0.2745, 1.0000];
    cyan   = [0.2706, 0.8784, 0.9686];
    green  = [0.4667, 0.6745, 0.1882];


    thmHom = struct;
    thmHom.special  = {'S', 'HB', 'MX', 'BP', 'NSS_', 'NSF_', 'DRS_', 'DRU_', 'TLS_', 'TLR_', ...
                       'NDS_', 'NDU_', 'OFS_', 'OFU_', 'IFS_', 'IFU_', 'NCH_', 'ENDSt_', 'ENDUn_', 'RES_', 'EqType_'};
    thmHom.ustab = 'RES';
    thmHom.ustabfun = @(x) 1 + (~isnan(x) & x > 0) + 2 * isnan(x);
    % thmHom.usept = thmHom.special;
    
    % Define marker styles for different homoclinic points
    thmHom.NSF_ = {'diamond', 'MarkerFaceColor', cyan, 'MarkerSize', 8};    
    thmHom.RES_ = {'diamond', 'MarkerFaceColor', purple, 'MarkerSize', 8};
    thmHom.EqType_ = {'diamond', 'MarkerFaceColor', purple, 'MarkerSize', 8};
    thmHom.NSS_ = {'diamond', 'MarkerFaceColor', purple, 'MarkerSize', 8};
    thmHom.BT_  = {'square', 'MarkerFaceColor', 'black', 'MarkerSize', 8};
    thmHom.S    = {'*', 'MarkerFaceColor', 'black', 'MarkerSize', 8};
    thmHom.HB   = {'*', 'MarkerFaceColor', 'black', 'MarkerSize', 8};
    thmHom.DRU_ = {'x', 'MarkerFaceColor', red, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.DRS_ = {'x', 'MarkerFaceColor', blue, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.NDS_ = {'+', 'MarkerFaceColor', red, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.NDU_ = {'+', 'MarkerFaceColor', blue, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.TLS_ = {'o', 'MarkerFaceColor', red, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.TLR_ = {'o', 'MarkerFaceColor', blue, 'MarkerSize', 10, 'linewidth', 2};
    thmHom.IFU_ = {'pentagram', 'MarkerFaceColor', red, 'MarkerSize', 8};
    thmHom.IFS_ = {'pentagram', 'MarkerFaceColor', blue, 'MarkerSize', 8};
    thmHom.OFU_ = {'hexagram', 'MarkerFaceColor', red, 'MarkerSize', 8};
    thmHom.OFS_ = {'hexagram', 'MarkerFaceColor', blue, 'MarkerSize', 8};
    thmHom.MX = {'o', 'MarkerEdgeColor', 'black', 'MarkerSize', 8};
    thmHom.BP = {'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 8};
    thmHom.lspec = {{'-', 'LineWidth', 1.5, 'color', green}, {'b-.', 'LineWidth', 1.5, 'color', green}};
end
