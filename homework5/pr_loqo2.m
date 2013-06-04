function [x,y] = pr_loqo2(c, H, A, b, l, u)
%[X,Y] = PR_LOQO2(c, H, A, b, l, u)
%
%loqo solves the quadratic programming problem
%
%minimize   c' * x + 1/2 x' * H * x
%subject to A'*x = b
%           l <= x <= u
%
% Dimensions: c : N-column vector
%             H : NxN matrix
%             A : N-row vector
%             b : real number
%             l : N-column vector
%             u : N-column vector
% 
%             x : N-column vector
%             y : Objective value
%             
%for a documentation see R. Vanderbei, LOQO: an Interior Point Code
%                        for Quadratic Programming
margin = 0.05; 
bound  = 100; 
sigfig_max = 8; 
counter_max = 50;
[m, n] = size(A); 
H_x    = H; 
H_diag = diag(H);
b_plus_1 = 1; 
c_plus_1 = norm(c) + 1;
one_x = -ones(n,1); 
one_y = -ones(m,1);

for i = 1:n 
    H_x(i,i) = H_diag(i) + 1; 
end;

H_y = eye(m); 
c_x = c; 
c_y = 0;
R = chol(H_x); 
H_Ac = R \ ([A; c_x'] / R)';

H_A = H_Ac(:,1:m); 
H_c = H_Ac(:,(m+1):(m+1));

A_H_A = A * H_A; A_H_c = A * H_c;
H_y_tmp = (A_H_A + H_y); y = H_y_tmp \ (c_y + A_H_c);
x = H_A * y - H_c; g = max(abs(x - l), bound);
z = max(abs(x), bound); t = max(abs(u - x), bound);
s = max(abs(x), bound); mu = (z' * g + s' * t)/(2 * n);
sigfig = 0; counter = 0; alfa = 1;
while ((sigfig < sigfig_max) * (counter < counter_max)),
  counter = counter + 1; H_dot_x = H * x;
  rho = - A * x + b; nu = l - x + g; tau = u - x - t;
  sigma = c - A' * y - z + s + H_dot_x;
  gamma_z = - z; gamma_s = - s;
  x_dot_H_dot_x = x' * H_dot_x;
  primal_infeasibility = norm([tau; nu]) / b_plus_1;
  dual_infeasibility = norm([sigma]) / c_plus_1;
  primal_obj = c' * x + 0.5 * x_dot_H_dot_x;
  dual_obj = - 0.5 * x_dot_H_dot_x + l' * z - u' * s + b'*y; %%%
  old_sigfig = sigfig;
  sigfig = max(-log10(abs(primal_obj - dual_obj)/(abs(primal_obj) + 1)), 0);
  hat_nu = nu + g .* gamma_z ./ z; hat_tau = tau - t .* gamma_s ./ s;
  d = z ./ g + s ./ t;
  for i = 1:n H_x(i,i) = H_diag(i) + d(i); end;
  H_y = 0;  c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
  c_y = rho; R = chol(H_x); H_Ac = R \ ([A; c_x'] / R)';
  H_A = H_Ac(:,1:m); H_c = H_Ac(:,(m+1):(m+1));
  A_H_A = A * H_A; A_H_c = A * H_c; H_y_tmp = (A_H_A + H_y);
  delta_y = H_y_tmp \ (c_y + A_H_c); delta_x = H_A * delta_y - H_c;
  delta_s = s .* (delta_x - hat_tau) ./ t;
  delta_z = z .* (hat_nu - delta_x) ./ g;
  delta_g = g .* (gamma_z - delta_z) ./ z;
  delta_t = t .* (gamma_s - delta_s) ./ s;
  gamma_z = mu ./ g - z - delta_z .* delta_g ./ g;
  gamma_s = mu ./ t - s - delta_s .* delta_t ./ t;
  hat_nu = nu + g .* gamma_z ./ z;
  hat_tau = tau - t .* gamma_s ./ s;
  c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
  c_y = rho; H_Ac = R \ ([A; c_x'] / R)';
  H_A = H_Ac(:,1:m); H_c = H_Ac(:,(m+1):(m+1));
  A_H_A = A * H_A; A_H_c = A * H_c;
  H_y_tmp = (A_H_A + H_y); delta_y = H_y_tmp \ (c_y + A_H_c);
  delta_x = H_A * delta_y - H_c; delta_s = s .* (delta_x - hat_tau) ./ t;
  delta_z = z .* (hat_nu - delta_x) ./ g;
  delta_g = g .* (gamma_z - delta_z) ./ z;
  delta_t = t .* (gamma_s - delta_s) ./ s;
  alfa = - 0.95 / min([delta_g ./ g; delta_t ./ t;
                      delta_z ./ z; delta_s ./ s; -1]);
  mu = (z' * g + s' * t)/(2 * n);
  mu = mu * ((alfa - 1) / (alfa + 10))^2;
  x = x + delta_x * alfa; g = g + delta_g * alfa;
  t = t + delta_t * alfa; y = y + delta_y * alfa;
  z = z + delta_z * alfa; s = s + delta_s * alfa;
end