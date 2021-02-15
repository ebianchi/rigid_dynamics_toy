function [D] = D_rod(q, l)
%D_rod : Generates the (n-by-p*k) D matrix for the rod toy problem,
%        evaluated at the provided value of q.
%
% The D matrix is n-by-p*k, where the first n-by-k chunk corresponds to the
% k configuration space vectors (n-by-1) that point in the tangential
% directions to the first point of contact.  The pth n-by-k chunk
% corresponds to the k configuration space vectors (n-by-1) that point in
% the tangential directions to the pth point of contact.

% The following block of code can be uncommented in order to show how the
% way to evaluate the matrix D is done.

%     syms l m r I mu g x y theta dx dy dtheta real;
% 
%     % two points of contact
%     q = [x; y; theta];
%     dq = [dx; dy; dtheta];
% 
%     % tangent vector to contact in point space (x,y)
%     t_hat1 = [1; 0];
%     t_hat2 = -t_hat1;
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
%     tq_c1_1 = t_hat1' * J_dq_to_dp1
%     tq_c1_2 = t_hat2' * J_dq_to_dp1
%     tq_c2_1 = t_hat1' * J_dq_to_dp2
%     tq_c2_2 = t_hat2' * J_dq_to_dp2
%
%     % obtain a symbolic expression for D as a function of elements in q
%     D = [tq_c1_1', tq_c1_2', tq_c2_1', tq_c2_2']

    x = q(1);
    y = q(2);
    theta = q(3);
    
    D = [1, -1, 1, -1;
         0,  0, 0,  0;
         -l/2*sin(theta), l/2*sin(theta), l/2*sin(theta), -l/2*sin(theta)];
end

