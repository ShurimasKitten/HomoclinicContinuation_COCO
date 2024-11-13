function [data, y] = homoclinicTestFunction(prob, data, u)
    % HIGHPO_CHECKTESTFUNC Test functions for homoclinic bifurcations in 3D systems
    %
    % This function evaluates various test conditions related to homoclinic
    % bifurcations in a three-dimensional dynamical system. It assumes a 3D
    % solution space and does not compute the three leading eigenvalues or
    % the Bi-focus condition.
    %
    % ARGS:
    %   - prob: Standard COCO (Continuation COntinuation) problem structure.
    %   - data: Data structure containing problem-specific information.
    %   - u: Vector containing state variables, parameters, and time interval.
    %
    % RETURNS:
    %   - data: Updated data structure (if modified within the function).
    %   - y: Vector of evaluated test functions.
    %
    % TODO:
    %   - Implement manifold conditions.
    %   - Compute orientability index.
    %
    % BUG REPORT:
    %  
    
    %% Initialize Constants and Extract Variables
    % Machine error tolerance to handle numerical precision issues
    ETOL = 1e-10; 
    ETOL_M = 1e-6;  % Lower tolerance for manifold-related conditions
    
    % Extract state and parameter dimensions from the data structure
    % TODO: Improve this to dynamically determine dimensions instead of assuming
    xdim = data.xdim;
    pdim = data.pdim;
    
    % Extract equilibrium and manifold points from the input vector 'u'
    x_ss = u(2*xdim+1 : 3*xdim);           % Equilibrium point (steady state)
    x_1 = u(1:xdim);                        % Final point on the stable manifold (x(1))
    x_0 = u(xdim+1:2*xdim);                 % Initial point on the unstable manifold (x(0))
    T = u(end);                              % Time interval length
    params = u(3*xdim+1 : 3*xdim+pdim);      % Parameter values
    f = data.f; % vec.f 
    
    % Use try/catch to handle potential errors gracefully, especially if
    % continuation does not stop 
    % try
        %% Eigenspace Setup
        % Retrieve eigenvectors and eigenvalues associated with the equilibrium
        [v_un, v_st, v_unT, v_stT, ev_un, ev_st, saddleVal] = getTestFunctionEvals(x_ss', params, f);
   
        % Extract the first (leading) stable and unstable eigenvalues
        firstSt = ev_st(end);   % Leading stable eigenvalue
        firstUn = ev_un(1);     % Leading unstable eigenvalue
        
        % Extract the first (transposed) stable and unstable eigenvectors
        % v_st_First = v_st(:, end);
        % v_un_First = v_un(:, 1);
        v_stT_First = v_stT(:, end);
        v_unT_First = v_unT(:, 1);
        
        % Determine the second stable and unstable eigenvalues if available
        if size(v_st, 2) >= 2 && size(v_un, 2) >= 2
            secondSt = ev_st(end-1);
            secondUn = ev_un(2);
        elseif size(v_un, 2) >= 2 && size(v_st, 2) < 2
            secondUn = ev_un(2);
            secondSt = NaN;
        elseif size(v_st, 2) >= 2 && size(v_un, 2) < 2
            secondSt = ev_st(end-1);
            secondUn = NaN;
        end
        
        %% Local Conditions
        %%% Check equilibrium type
        % Resonance condition based on saddle value
        ResTest = saddleVal;
        
        % Initialize equilibrium type (+/-)
        %  1  : Real leading eigenvalues
        %  2  : Saddle-focus with 2D stable manifold
        %  3  : Saddle-focus with 2D unstable manifold
        %  4  : Bi-focus
        EqType = 1; 

        % Determine the type of equilibrium based on eigenvalues
        % Change sign each time a change is detected to 'cross zero' 
        if xdim >= 2 && size(v_st, 2) >= 2 && ~isreal(secondSt) && ~isreal(firstSt) && isreal(secondUn) && isreal(firstUn)  
            % Saddle-focus with 2D stable manifold
            EqType = sign(EqType)*2;  
        elseif xdim >= 2 && size(v_un, 2) >= 2 && ~isreal(secondUn) && ~isreal(firstUn) && isreal(secondSt) && isreal(firstSt)
            % Saddle-focus with 2D unstable manifold
            EqType = sign(EqType)*3;  
        elseif xdim >= 4 && size(v_un, 2) >= 2 && size(v_st, 2) >= 2 && all(~isreal(ev_un(1:2))) && all(~isreal(ev_st(end-1:end)))      
            % Bi-focus
            EqType = sign(EqType)*4; 
        end
        
        %%% Resonance Conditions
        if abs(ResTest)<ETOL_M && abs(EqType)==1
            neutralReal = 0;  
            fprintf("Neutral saddle type with real eigenvalues. \n");
        elseif abs(EqType)==1
            neutralReal = ResTest; 
        else
            neutralReal = 1;
        end

        if abs(ResTest)<ETOL_M && abs(EqType)==2
            neutralSF = 0; 
            fprintf("Saddle-focus transition (stable) i.e., Belnikov case 1.\n");
        elseif abs(ResTest)<ETOL_M && abs(EqType)==3
            neutralSF = 0; 
            fprintf("Saddle-focus transition (unstable) i.e., Belnikov case 1.\n");
        elseif abs(EqType)==2 || abs(EqType)==3
            neutralSF = ResTest;
        else 
            neutralSF = 1;
        end
       
        %%% Double Real Leading Eigenvalues        
        % Double real leading stable eigenvalues
        if abs(imag(firstSt)) < ETOL_M
            doubleR_St = (real(firstSt) - real(secondSt))^2;
        else
            doubleR_St = -(imag(firstSt) - imag(secondSt))^2;
        end
        
        % Double real leading unstable eigenvalues
        if abs(imag(firstUn)) < ETOL_M
            doubleR_Un = (real(firstUn) - real(secondUn))^2;
        else
            doubleR_Un = -(imag(firstUn) - imag(secondUn))^2;
        end
        
        % Neutrally Divergent Saddle-Focus (Stable)        
        NDSF_st = real(firstSt) + real(secondSt) + real(firstUn);
        
        % Neutrally Divergent Saddle-Focus (Unstable)        
        NDSF_un = real(firstSt) + real(secondUn) + real(firstUn);
        
        %%% Three Leading Eigenvalues. This encomposes the 3DLS senario. 
        % Unstable
        if length(ev_un)>=3
            threeLeadingUn = real(ev_un(1)) + real(ev_un(3));
        else
            threeLeadingUn = 0;
        end
        
        % Stable
        if length(ev_st)>=3 
            threeLeadingSt = real(ev_st(end)) - real(ev_st(end-2));
        else
            threeLeadingSt = 0;
        end 

        %%% Possible end of curve
        END = real(firstSt)*real(secondSt);
        
        %% Orbit Conditions
        
        %%% Stable flips
        if ~isreal(firstSt)
            flipS = abs(exp(-real(firstSt)*T) * dot(real(v_stT_First), x_1 - x_ss))^2 + abs(exp(-real(firstSt)*T) * dot(imag(v_stT_First), x_1 - x_ss))^2;
        elseif isreal(firstSt) 
            flipS = abs(exp(-real(firstSt)*T) * dot(real(v_stT_First), x_1 - x_ss));
        end
            
        % This handles the discontinious nature of the test functino
        if flipS < ETOL_M
            flipS = 0;
        end
        
        %%% Unstable flips
        if ~isreal(firstUn)
             flipU = abs(exp(real(firstUn)*T) * dot(real(v_unT_First), x_0 - x_ss))^2 + abs(exp(real(firstUn)*T) * dot(imag(v_unT_First), x_0 - x_ss))^2;
        elseif isreal(firstUn) 
             flipU = abs(exp(real(firstUn)*T) * dot(real(v_unT_First), x_0 - x_ss));
        end

        % This handles the discontinious nature of the test functino
        if flipU < ETOL_M
            flipU = 0;
        end
        
        %%% Compute Coefficients for Orbit Flip Detection        
        if flipS==0
            % Compute coefficients alpha and beta based on stable manifold conditions
            vss     = real(ev_st(1));
            vsLead  = ev_st(end);
            vuLead  = ev_un(1);
            alpha   = -vss / vuLead;
            beta    = -vsLead / vuLead;
        elseif flipU==0
            % Compute coefficients alpha and beta based on unstable manifold conditions
            vuu     = real(ev_un(end));
            vuLead  = ev_un(1);
            vsLead  = ev_st(end);
            alpha   = -vuu / vsLead;
            beta    = -vuLead / vsLead;
        end
        
        %%% Print Orbit Flip Type
        if flipU==0 || flipS==0
            if beta > 1
                fprintf("ORBIT FLIP TYPE A: BELNIKOV => One-sided N-homoclinics!\n");
            elseif beta < 1 && alpha > 1
                fprintf("ORBIT FLIP TYPE B: Homoclinic-doubling! => 1-PD, 2-PO, 2-Hom!\n");
            elseif alpha < 1
                fprintf("ORBIT FLIP TYPE C: Scary case! \n");
            else
                fprintf("ERROR IN ORBIT FLIP TYPE CLASSIFICATION!\n")
            end
        end


        %% Manifold Conditions (TODO)        
        % Placeholder variables for manifold conditions
        inlMS = -1;
        inlMU = -1;

        %% Output Test Functions
        y = [neutralReal, neutralSF, doubleR_St, doubleR_Un, threeLeadingSt, threeLeadingUn, ...
             NDSF_st, NDSF_un, flipS, flipU, inlMS, inlMU, ResTest, EqType, END]';
        
    % catch 
    %     %% Something went wrong!
    %     % If an error occurs during computation, set all test functions to NaN
    %     % This may indicate the end of a continuation branch or other issues
    % 
    %     % In-place error handling!
    %     if data.counter > 1
    %         error("Something went wrong! Error in selecting EVALS! Possible end of run?")
    %     end
    %     fprintf("Something went wrong! Error in selecting EVALS! Possible end of run?\n\n")
    %     y = NaN * ones(1, 14);
    %     data.counter = data.counter + 1;
    % end
end
