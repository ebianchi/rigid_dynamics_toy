function [G] = G_rod(q, M)
%G_rod : Generates the (n-by-1) G vector for the rod toy problem, evaluated
%        at the provided value of q.

% The following block of code can be uncommented in order to show how the
% way to evaluate the vector G is done.
%
%     syms m g x y theta real;
%     q = [x; y; theta];
% 
%     % get the potential energy at q
%     V = m*g*[0, 1, 0]*q;  % select the y component of q
%     
%     % convert the potential energy to gravity vector
%     G = jacobian(V, q)

    g = 9.81;   % [m/s^2]
    m = M(1,1);  % [kg]
    
    G = [0; m*g; 0];
end

