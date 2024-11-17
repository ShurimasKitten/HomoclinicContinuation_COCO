function [v_un, v_st] = buildExplBC(x0, p0, data)
    % Finds the stable and unstable eigenvectors of the Jacobian matrix for the
    % x0_neg non-trivial equilibrium point.
    %

    global f
    
    J = f.dfdx(x0, p0);
    % Calculate eigenvalues
    [eigvec, eigval] = eig(J);

    % Sort
    [eigval,ind] = sort(diag(eigval));
    eigvec = eigvec(:,ind);
    
    % Inidices for unstable eigenvector (eigval > 0)
    unstable_index = find(eigval > 0);
    % Indices for stable eigenvectors (eigval < 0)
    stable_index = find(eigval < 0);

    % (un)stalb eigenvalues
    ev_un = eigval(unstable_index,:);
    ev_st = eigval(stable_index,:); 

    % (un)stable eigenvectors
    v_un = eigvec(:, unstable_index);
    v_st = eigvec(:, stable_index);

    % Normalise e-vec
    for i=1:length(unstable_index)
      v_un(:,i) = v_un(:,i)/norm(v_un(:,i));
    end

    for i=1:length(stable_index)
      v_st(:,i) = v_st(:,i)/norm(v_st(:,i));
    end
    v_un = v_un';
    v_st = v_st';
end

