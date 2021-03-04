clear
syms x y;

% Cost function and its derivatives
f(x, y) = (1 - x)^2  + 100*(y - x^2)^2;
g = gradient(f);
A = hessian(f);

% Convert cost funciton its derivatives to numerical functions
f = matlabFunction(f);
g = matlabFunction(g);
A = matlabFunction(A);

maxIts = 200;       % Maximum iterations allowed
epsilon = 10^(-8);  % Stopping tolerance
B = [2 0; 0 0.5];   % B defines the shape of the elliptical trust region
Delta = 3;          % Initial "radius" of trust region
x = [0; 0];         % Initial solution guess

% M_tilde pencil defintion
ML = @(x,y) [-B, A(x,y); A(x,y), -g(x,y)*(g(x,y)')/Delta^2];  % LHS of M-tilde
MR = [zeros(2,2), -B; -B, zeros(2,2)];                        % RHS of M-tilde

% Algorithm 5.1 implementation
for i = 1:maxIts
    if (norm(g(x(1), x(2))) < eps)
        break;
    end
    
    % Compute Newton step
    p0 = -A(x(1),x(2)) \ g(x(1),x(2));
    
    % Ignore Newton step if it is outside of the trust region
    if ~(BNorm(p0, B) < Delta)
        p0 = [0;0];
    end
    
    % Find the rightmost eigenvalue/vector of the pencil M_tilde
    [y, lambda] = eigs(ML(x(1),x(2)), MR, 1, 'largestreal');
    y = real(y);
    y1 = y(1:2);
    y2 = y(3:4);
    
    % Compute boundary step from rightmost eigenvector
    p1 = -sign(g(x(1),x(2))' * y2) * Delta * y1 / BNorm(y1, B);
    
    % Choose the best step between Newton and boundary steps
    if (f(x(1) + p0(1), x(2) + p0(2)) < f(x(1) + p1(1), x(2) + p1(2)))
        x = x + p0;
    else
        x = x + p1;
    end
end

x       % Minimizer of f


% Elliptical B-Norm definition
function val = BNorm(x, B)
    val = sqrt(x'*B*x);
end