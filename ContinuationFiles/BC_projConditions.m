function [data, y] = BC_projConditions(prob, data, u)
    % BC_projConditions defines boundary conditions for homoclinic continuation.
    %
    % This function sets up boundary conditions at the start (U(0)) and end (U(1))
    % of each continuation segment to ensure that they lie along the tangent
    % spaces of the left-hand side (LHS) and right-hand side (RHS) equilibria.
    %
    % This function is called at every continuation step.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - data: Data structure containing necessary information.
    %   - u: State and parameter variables concatenated into a single vector.
    %
    % Outputs:
    %   - data: Updated data structure (unchanged in this function).
    %   - y: Boundary condition residuals enforcing the projection conditions.
    
    % Extract state- and parameter-space dimensions from data
    xdim = data.xdim;  % Dimension of the state space
    pdim = data.pdim;  % Dimension of the parameter space

    % Parse the input vector 'u' into individual components
    x_1 = u(1:xdim);                     % Point in the stable manifold (Ws)
    x_0 = u(xdim+1:2*xdim);              % Point in the unstable manifold (Wu)
    x_ss = u(2*xdim+1 : 3*xdim);         % Equilibrium point
    params = u(3*xdim+1 : 3*xdim+pdim);  % System parameters

    % Obtain projection matrices based on the equilibrium and parameters
    [Lu, Ls] = computeManifoldProjectionMatrices(x_ss, params, data);

    % Compute boundary conditions for the unstable manifold
    unstableBC = Lu * (x_0 - x_ss);

    % Compute boundary conditions for the stable manifold
    stableBC   = Ls * (x_1 - x_ss);

    % Combine stable and unstable boundary conditions into output
    y = [stableBC; unstableBC];
end
