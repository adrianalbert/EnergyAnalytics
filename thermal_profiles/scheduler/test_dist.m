a = [1 2 3 4]';
b = [3 4 5 2]';
n = ones(4,1);
B = 3;
p = 4;

cvx_begin
    variable u(p);
    minimize(norm(a - u .* b,2));
    subject to
        u' * n <= B;
        u >= 0;
        u <= 1;
cvx_end