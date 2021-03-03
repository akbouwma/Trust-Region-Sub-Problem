clear
syms x y;

% We match the papers naming conventions
f(x, y) = (1 - x)^2  + 100*(y - x^2)^2;

B = eye(2);  % In the paper, B is the scaling matrix, we probably want to change this from the identity to something else
Delta = 0.1;

% find the gradient and hessian of f
g = gradient(f);
A = hessian(f);  % In the paper, A is the Hessian

x0 = [1, 1];



