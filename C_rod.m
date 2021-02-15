function [C] = C_rod(q, v, M)
%C_rod : Generates the (n-by-n) C matrix for the rod toy problem, evaluated
%        at the provided values of q and v.

% The following block of code can be uncommented in order to show how the
% way to evaluate the matrix C is done.
%
%     syms l m r I mu g x y theta dx dy dtheta real;
% 
%     M = diag([m, m, I]);
%     v = [dx; dy; dtheta];
%     q = [ x;  y;  theta];
% 
%     mv_grad_q = jacobian(M*v, q);
%     C = mv_grad_q - 1/2 * mv_grad_q'
    
    C = [0, 0, 0; 0, 0, 0; 0, 0, 0];
end

