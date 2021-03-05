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
B = [1,1; 1,1];     % B defines the shape of the elliptical trust region
B = B*B';           % Make sure B is SPD
Delta = 1;          % Initial "radius" of trust region
x = [0; 0];         % Initial solution guess

% M_tilde pencil defintion
M0 = @(x,y) [-B, A(x,y); A(x,y), -g(x,y)*(g(x,y)')/Delta^2]; % LHS of M-tilde
M1 = -[zeros(2,2), B; B, zeros(2,2)];                        % RHS of M-tilde

% Algorithm 5.1 implementation
for i = 1:maxIts
    if (norm(g(x(1), x(2))) < eps)
        break;
    end
    
    % Compute Newton step
    p0 = -A(x(1),x(2)) \ g(x(1),x(2));
    
    % Ignore Newton step if it is outside of the trust region
    if ~(BNorm(p0, B) < Delta)
        p0TEST = false;        % Ignore Newton step, outside of region
    else
        p0TEST = true;         % Allow Newton step this iteration
    end
    
    % Find the rightmost eigenvalue/vector of the pencil M_tilde
    [y, lambda] = eigs(M0(x(1),x(2)), M1, 1, 'largestreal');
    p1 = y(3:4);
    p1 = p1/BNorm(p1,B)*Delta;
    
    % Choose the best step between Newton and boundary steps
    if (p0TEST && f(x(1) + p0(1), x(2) + p0(2)) < f(x(1) + p1(1), x(2) + p1(2)))
        x = x + p0;
    else
        x = x + p1;
        p0TEST = true;
    end
end

i       % Iterations used
x       % Minimizer of f


% Elliptical B-Norm definition
function val = BNorm(x, B)
    val = sqrt(x'*B*x);
end