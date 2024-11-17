function [hom_data, y] = BC_phase(prob, hom_data, u)
    % BC_hom_phase defines the discretized phase condition for 'coll' toolbox continuation.
    %
    
    % Compute the phase condition
    y = hom_data.xp0 * u;
end
