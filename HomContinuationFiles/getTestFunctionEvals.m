function [v_un, v_st, v_unT, v_stT, ev_un, ev_st, sadVal] = getTestFunctionEvals(x0, p0, f)
    % Finds the stable and unstable eigenvectors of the Jacobian matrix for the
    % x0_neg non-trivial equilibrium point.
    %
    
    function [v_un, v_st, ev_un, ev_st] = findEspaces(J)
        % Calculate eigenvalues
        [eigvec, eigval] = eig(J);

        % Sort
        eigval = diag(eigval); % Convert to a vector if needed
        [~, ind] = sort(real(eigval)); % Sort by the real part of the eigenvalues
        eigval = eigval(ind);  % Rearrange eigenvalues
        eigvec = eigvec(:,ind);
            
        % Inidices for unstable eigenvector (eigval > 0)
        unstable_index = eigval > 0;
        % Indices for stable eigenvectors (eigval < 0)
        stable_index = eigval < 0;
    
        % (un)stalb eigenvalues
        ev_un = eigval(unstable_index,:);
        ev_st = eigval(stable_index,:); 
       
        % (un)stable eigenvectors
        v_un = eigvec(:, unstable_index);
        v_st = eigvec(:, stable_index);
        
        % Normalise 
        for i = 1:size(v_un,2)
            v_un(:,i) = v_un(:,i) / norm(v_un(:,i));
        end
        for i = 1:size(v_st,2)
            v_st(:,i) = v_st(:,i) / norm(v_st(:,i));
        end
    end

    % Espaces of J
    J = f.dfdx(x0', p0);
    [v_un, v_st, ev_un, ev_st] = findEspaces(J);

    % Espaces of J^T
    JT = f.dfdx(x0', p0)';
    [v_unT, v_stT, ~, ~] = findEspaces(JT);
    
    % Compute saddle-val -- leading unstable over leading stable
    sadVal = real(ev_un(1)) + real(ev_st(end));
end