function [N] = N_rod(q, v, l)
%N_rod : Generates the (n-by-p) N matrix for the rod toy problem, evaluated
%        at the provided value of q.
%
% The N matrix is n-by-p, where each n-by-1 column corresponds to the
% vector normal of the pth contact.

% The following block of code can be uncommented in order to show how the
% way to evaluate the matrix N is done.

%     syms l m r I mu g x y theta dx dy dtheta real;
% 
%     % two points of contact
%     q = [x; y; theta];
%     dq = [dx; dy; dtheta];
% 
%     % normal vector to contact in point space (x,y)
%     nhat = [0; 1];
% 
%     % positions of contact points in point space (x,y)
%     p_c1  = [x + l/2*cos(theta); y + l/2*sin(theta) - r];
%     p_c2  = [x - l/2*cos(theta); y - l/2*sin(theta) - r];
% 
%     % velocities of contact points in point space (x,y)
%     dp_c1 = [dx - l/2*sin(theta)*dtheta; dy + l/2*cos(theta)*dtheta];
%     dp_c2 = [dx + l/2*sin(theta)*dtheta; dy - l/2*cos(theta)*dtheta];
% 
%     % get the jacobians to convert point velocities to configuration
%     % velocities
%     J_dq_to_dp1 = jacobian(dp_c1, dq);
%     J_dq_to_dp2 = jacobian(dp_c2, dq);
% 
%     % get the contact normal in configuration space (x,y,theta)
%     nq_c1 = nhat' * J_dq_to_dp1;
%     nq_c2 = nhat' * J_dq_to_dp2;
% 
%     % obtain a symbolic expression for N as a function of elements in q
%     n_of_q = [nq_c1', nq_c2']

    x = q(1);
    y = q(2);
    theta = q(3);

    N = [0, 0;
         1, 1;
         l/2 * cos(theta), -l/2 * cos(theta)];
end

