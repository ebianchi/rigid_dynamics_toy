function [phis] = phis_rod(q, l, r)
%phiss_rod : Generates the (p-by-1) phis vector for the rod toy problem,
%            evaluated at the provided value of q.

% The first contact point is the right one, and the second is the left one.

    x = q(1);
    y = q(2);
    theta = q(3);

    phis = [y + l/2*sin(theta) - r;
            y - l/2*sin(theta) - r];
end

