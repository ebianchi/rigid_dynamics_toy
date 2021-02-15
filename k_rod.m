function [k] = k_rod(q, v, M)
%k_rod : Generates the (n-by-1) k vector for the rod toy problem, evaluated
%        at the provided value of q and v.
    
    k = -C_rod(q,v,M)*v - G_rod(q, M);
end

