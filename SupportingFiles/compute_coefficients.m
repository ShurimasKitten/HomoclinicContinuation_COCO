
function c = compute_coefficients(x_ss, Ls, x1)
    % COMPUTE_COEFFICIENTS Computes the coefficients for expressing x1 in terms of x_ss and Ls
    %
    %   c = COMPUTE_COEFFICIENTS(x_ss, Ls, x1) returns the coefficients vector c
    %   such that:
    %       x1 = x_ss + Ls * c
    %
    %   Inputs:
    %       x_ss - (n x 1) equilibrium point vector
    %       Ls   - (n x k) matrix of stable eigenvectors (each column is an eigenvector)
    %       x1   - (n x 1) point to be expressed as x_ss + Ls * c
    %
    %   Output:
    %       c    - (k x 1) coefficients vector

    % Compute the difference vector
    x_diff = x1 - x_ss;
       
    c = Ls \ x_diff;
    
    % Alternatively, use the pseudo-inverse for more stability in some cases
    % c = pinv(Ls) * x_diff;
    c = c/norm(c);

    % Optional: Verify the solution (uncomment if needed)
    % x_reconstructed = x_ss + Ls * c;
    % residual = norm(x1 - x_reconstructed);
    % fprintf('Residual: %e\n', residual);
end
