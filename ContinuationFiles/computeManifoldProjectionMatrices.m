function [Lu, Ls] = computeManifoldProjectionMatrices(x0, p0, data)
    % buildProjBC constructs projection matrices for stable and unstable manifolds.
    %
    % This function computes the stable (Ls) and unstable (Lu) projection matrices
    % based on the Jacobian matrix at the equilibrium point (x0) with parameters (p0).
    %
    % Inputs:
    %   - x0: Equilibrium point (state vector) at which the Jacobian is evaluated.
    %   - p0: Parameter vector corresponding to the equilibrium point.
    %   - data: Data structure containing necessary information, including eigendata.
    %
    % Outputs:
    %   - Lu: Projection matrix for the unstable manifold.
    %   - Ls: Projection matrix for the stable manifold.
    
    % Extract the system's function structure from data
    f = data.f;
    
    % Compute the Jacobian matrix at the equilibrium point x0 with parameters p0
    J = f.dfdx(x0, p0);
 
    % Calculate eigenvalues and eigenvectors of the Jacobian
    [eigvec, eigval] = eig(J);
    
    % Sort eigenvalues in ascending order and rearrange eigenvectors accordingly
    [eigval, ind] = sort(diag(eigval));
    eigvec = eigvec(:, ind);
    
    % Identify indices of unstable (eigval > 0) and stable (eigval < 0) eigenvalues
    unstable_index = eigval > 0;
    stable_index = eigval < 0;
    
    % Extract unstable and stable eigenvectors based on identified indices
    v_un = eigvec(:, unstable_index);
    v_st = eigvec(:, stable_index);
    
    % Normalize the columns of unstable and stable eigenvectors
    v_un = v_un ./ vecnorm(v_un, 2, 1); 
    v_st = v_st ./ vecnorm(v_st, 2, 1);
    
    % Compute the orthogonal complements of the unstable and stable eigenvectors
    v_un_star = null(v_un')';  % Orthogonal complement to unstable eigenvectors
    v_st_star = null(v_st')';  % Orthogonal complement to stable eigenvectors
    
    % Check if 'data.eigdata' is provided; if not, use the current orthogonal complements
    if nargin < 3 || isempty(data.eigdata)
        % Initial setup: assign orthogonal complements directly
        Lu = v_un_star;
        Ls = v_st_star;
    else
        % Use previously stored orthogonal complements from 'data.eigdata'
        V_ST_tilda = data.eigdata{1};  % Previously computed orthogonal complements for stable manifold
        W_UN_tilda = data.eigdata{2};  % Previously computed orthogonal complements for unstable manifold

        % Current orthogonal complements
        V_ST_star = v_st_star;
        W_UN_star = v_un_star;

        % Compute projection matrices using the method from Byen 1990
        U_s = (V_ST_tilda * V_ST_tilda') / (V_ST_star * V_ST_tilda');
        U_u = (W_UN_tilda * W_UN_tilda') / (W_UN_star * W_UN_tilda');
        
        % Construct the stable and unstable projection matrices
        Ls = U_s * V_ST_star;
        Lu = U_u * W_UN_star;
    end
end
