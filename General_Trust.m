%{
% Solve N-D Elliptical TRS
% Authors: Alan Bouwman, Caleb Jacobs
%}

% Clean workspace
clear
close

% Number of dimensions
n = 2;

% Cost function definitions
% Cost function
f = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2; % Cost function

% Gradient of cost function
G = @(x) [-2*(1 - x(1)) - 400 * x(1) * (x(2) - x(1)^2); ...
          200*(x(2) - x(1)^2)]; 
      
% Hessian of cost function
H = @(x) [2 + 800*x(1)^2 - 400*(x(2) - x(1)^2), -400*x(1); ...
-400*x(1),	200];


% Trust region parameters
maxIts  = 200;                   % Max iterations allowed for convergence
epsilon = 10^(-5);               % Stopping tolerance
Delta   = 1;                     % Initial trust region "radius"
scale   = 0.1;                   % Delta changer
x       = [100;5];                 % Initial minimizer guess

% Elliptical trust region matrix
B = getEllipticalMatrix(2, [1;1], 3, [1;-1]);

% Begin solving the TRS
for i = 1:maxIts
    % Compute current gradient
    g = G(x);
    
    % Check stopping criteria
    if (norm(g) < epsilon)
        break;
    end
    
    % Compute current hessian
    A = H(x);
    
    % Compute Newton step
    p0 = -A \ g;
    
    % Assume we will use Newton step
    useP0 = false;
    
    % See if Newton step is feasible (i.e. ||p0||_B < Delta)
    if BNorm(p0, B) < Delta
        useP0 = true;   % Newton step is feasible, so use it    
    end
    
    % Compute left matrix of M
    M0 = [-B A;A, -g*g'/Delta^2];
          
    % Compute right matrix of M
    M1 = [zeros(n) B;B zeros(n)];
    
    % Find the rightmost eigenvalue/vector of the pencil M_tilde
    [y, lambda] = eigs(M0, -M1, 1, 'largestreal');
    y1 = y(1:n);
    y2 = y(n+1:end);
    
    % Compute boundary step from rightmost eigenvector
    p1 = -sign(g' * y2) * Delta * y1 / BNorm(y1, B);
    
    % Compare Newton step to TR-boundary step to get optimal solution
    if (useP0 && f(x + p0) < f(x + p1))
        x = x + p0;
%         Delta = Delta * (1 - scale);
        fprintf("p0: %d, x = (%.1f,%.1f), Delta = %.2f\n", i, x(1), x(2), Delta)
    else
        x = x + p1;
%         Delta = Delta * (1 + scale);
        fprintf("p1: %d, x = (%.1f,%.1f), Delta = %.2f\n", i, x(1), x(2), Delta)
    end
end

B
x
i

% Elliptical norm
function a = BNorm(p, B)
    a = sqrt(p'*B*p);
end

% Generate elliptical region matrix
function B = getEllipticalMatrix(a, va, b, vb)
    va = va / norm(va);     % Normalize major axis vector
    vb = vb / norm(vb);     % Normalize minor axis vector
    
    P = [va vb];
    D = diag([1/a,1/b]).^2;
    
    B = P * D * P^(-1);
end