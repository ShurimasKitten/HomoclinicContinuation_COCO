function [hom_data, y] = BC_projConditions(prob, hom_data, u)
    % BC_projConditions defines boundary conditions 'y' for homoclinic continuation.
    %

    % Dim info
    xdim = hom_data.xdim; 
    pdim = hom_data.pdim;  
    
    % Parse the input vector 'u' into individual components
    x_1 = u(1:xdim);                     % Point in the stable manifold (Ws)
    x_0 = u(xdim+1:2*xdim);              % Point in the unstable manifold (Wu)
    x_ss = u(2*xdim+1 : 3*xdim);         % Equilibrium point
    params = u(3*xdim+1 : 3*xdim+pdim);  % System parameters

    % Obtain projection matrices based on the equilibrium and parameters
    [Lu, Ls] = computeManifoldProjectionMatrices(x_ss, params, hom_data);

    % Compute boundary conditions for the unstable manifold
    unstableBC = Lu * (x_0 - x_ss);

    % Compute boundary conditions for the stable manifold
    stableBC   = Ls * (x_1 - x_ss);

    % Combine stable and unstable boundary conditions into output
    y = [stableBC; unstableBC];
end
