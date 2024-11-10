function [data, y] = BC_phase(prob, data, u)
    % BC_hom_phase defines the discretized phase condition for 'coll' toolbox continuation.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - data: Data structure containing necessary information, including the
    %           phase condition vector `xp0`.
    %   - u: State and parameter variables concatenated into a single vector.
    %
    % Outputs:
    %   - data: Updated data structure (unchanged in this function).
    %   - y: Enforcing the discretized phase constraint.
    %
    
    % Compute the phase condition
    y = data.xp0 * u;
end
