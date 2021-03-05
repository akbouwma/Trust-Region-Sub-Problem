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
x       = [0;0];                 % Initial minimizer guess

% Elliptical trust region matrix
B = [1,2;3,4];
B = B*B'                         % Ensure B is SPD

% Newton step validation
useP0 = false;

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
    
    % See if Newton step is feasible (i.e. ||p0||_B < Delta)
    if BNorm(p0, B) < Delta
        useP0 = true;   % Newton step is feasible, so use it    
    end
    
    % Compute left matrix of M
    M0 = [-B, A; A, -g*g' / Delta^2];
          
    % Compute right matrix of M
    M1 = [zeros(n), B; B, zeros(n)];
    
    % Compute TR-boundary step
    [v, lam] = eigs(@(x)M0*x, 2*n, -M1, 1, 'lr');
    p1 = v(n+1:end);
    p1 = -p1 / BNorm(p1, B) * Delta;
    
    % Compare Newton step to TR-boundary step to get optimal solution
    if (useP0 && f(x + p0) < f(x + p1))
        x = x + p0;
    else
        x = x + p1;
    end
    useP0 = false;
end

x

% Elliptical norm
function a = BNorm(p, B)
    a = sqrt(p'*B*p);
end