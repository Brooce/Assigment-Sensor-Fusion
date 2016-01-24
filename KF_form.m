function [X, P, K] = KF_form(Xs_1, Xs_k, h0, P_r_filt_ratio, x_state_prev, P_cov_prev, F,G,Q,R)

%=KF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF);



S1 = Xs_1';                          % UAV starting position vector (x0; y0)
Sk = Xs_k';                          % Current UAV position vector (x; y)
alpha = P_r_filt_ratio;              % Alpha - power ratio (measurement) 

Xp = x_state_prev;                   % State from the previous iteration
Pp = P_cov_prev;





%% Extended Kalman Filter

%% Prediction

% Equation 1
X_s = F * Xp;


% Equation 2
P_s = F * Pp * F' + G * Q * G';       % Error covariance extrapolation


%% Correction

% nonlinear measurement equation:

x = X_s(1);
y = X_s(2);
sx1 = S1(1);
sy1 = S1(2);
sx = Sk(1);
sy = Sk(2);

num_h = (x - sx1)^2 + (y - sy1)^2 + h0^2;
den_h = (x - sx)^2 + (y - sy)^2 + h0^2;

h = num_h / den_h;

% Jacobian of nonlinear measurement equation:
H = [   ( 2*(x - sx1) * den_h - 2*(x - sx) * num_h ) / den_h^2      ( 2*(y - sy1) * den_h - 2*(y - sy) * num_h ) / den_h^2      ];


% Equation 3: Innovation
v = (alpha - h);


% Equation 4: Innovation
S = H * P_s * H' + R;


% Equation 5: Kalman gain
K = P_s * H' / S;




% Equation 6: State update
X = X_s + K * v;



% Equation 7: Error covariance update
P = (eye(2) - K * H) * P_s;


%pause(0.1)

clc;


%keyboard

