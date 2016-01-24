function [X, P, K] = KF_form(Xs_1, Xs_k, h0, P_r_filt_ratio, x_state_prev, P_cov_prev, F,G,Q,R)

%=KF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF);



S1 = Xs_1';                          % UAV starting position vector (x0; y0)
Sk = Xs_k';                          % Current UAV position vector (x; y)
alpha = P_r_filt_ratio;              % Alpha - power ratio (measurement) 

Xp = x_state_prev;                   % State from the previous iteration
Pp = P_cov_prev;





%% Unscented Kalman Filter


%% ================================= PREDICTION ======================================




%% Augmenting the matrices with process noise:

E_wk = [0;
        0];          % Assuming that the process noise is gaussian with zero mean?

%
Xpa = [Xp;
       E_wk];
   
Ppa = [      Pp         zeros(size(Pp,1))       ;
      zeros(size(Q,1))            Q           ];
  
Xpa_Pdim = Xpa(:, ones( size(Ppa,1), 1 ));


if size(Q) ~= size(Pp)
    error('Nierowne wymiary P i Q');
end

%}


% Deriving a set of (2L + 1) points:

L = numel(Xpa);





%% Weights

a = 10^-3;
k = 1;
b = 2;

l = a^2 * (L - k) - L;

Ws0 = l/(L+l);                          % It will always be negative!!! ???
Ws0 = 1 / (2 * L);                      % Let's try like this...
Wc0 = Ws0 + (1 - a^2 + b);
Ws = 1/2 * Ws0;
Wc = Ws;

%% Sigma points:

%keyboard

%X_psig = [   Xpa,    Xpa_Pdim + sqrt(L + l) * chol(Ppa)',    Xpa_Pdim - sqrt(L + l) * chol(Ppa)'  ];
X_psig = [   Xpa,    Xpa_Pdim + sqrt( (L + l) * Ppa ),    Xpa_Pdim - sqrt( (L + l) * Ppa )  ];




%% Propagating through the state transition function and combining again:


for ( i = 1 : (2*L+1) )

    % Sigma points are propagated through F:
    X_sig(:,i) = X_psig(:,i);                   % F doesn't change anything but dimensions don't agree

end




% The weighted sigma points are recombined:
X_s = Ws * X_sig;
X_s(:,1) = Ws0 * X_sig(:,1);
X_s = sum(X_s, 2);


% To produce the covariance:
for ( i = 1 : (2*L+1) )
    P_sig(:,:,i) = Wc * ( X_sig(:,i) - X_s )*( X_sig(:,i) - X_s )';
end

P_s(:,:,1) = Wc0 * ( X_sig(:,1) - X_s )*( X_sig(:,1) - X_s )';
P_s = sum(P_sig, 3);


%keyboard





%% =================================== UPDATE ========================================

%% Augmenting the matrices with process noise:

E_v = zeros(size(R,1), 1);
%keyboard


Xsa = [X_s; E_v];
Psa = [      P_s         zeros(size(P_s, 1), 1)       ;
      zeros(1, size(P_s, 2))            R             ];
  
Xsa_Pdim = Xsa(:, ones( size(Psa,1), 1 ));







L = numel(Xsa);
l = a^2 * (L - k) - L;
Ws0 = l/(L+l);
Ws0 = 1 / (2 * L);
Wc0 = Ws0 + (1 - a^2 + b);
Ws = 1/2 * Ws0;
Wc = Ws;


%keyboard;
X_sig = [   Xsa,    Xsa_Pdim + sqrt( (L + l)*Psa ),    Xsa_Pdim - sqrt( (L + l)*Psa )  ];
%X_sig = [   Xsa,    Xsa_Pdim + sqrt(L + l) * chol(Psa)',    Xsa_Pdim - sqrt(L + l) * chol(Psa)'  ];

%keyboard


%% Propagating through the measurement function:


sx1 = S1(1);
sy1 = S1(2);

sx = Sk(1);
sy = Sk(2);

    
for ( i = 1 : (2*L+1) )
    
    
    x = X_sig(1, i);
    y = X_sig(2, i);

    num_h = (x - sx1)^2 + (y - sy1)^2 + h0^2;
    den_h = (x - sx)^2 + (y - sy)^2 + h0^2;

    gamma(i) = num_h / den_h;
    
end




%% Recombining to produce the predicted measurement and covariance:

z = Ws * gamma;
z(:,1) = Ws0 * gamma(:,1);
z = sum(z, 2);
%keyboard

for (i = 1 : (2*L+1) )
    %keyboard
    P_z(:,:,i) = Wc * ( gamma(i) - z )*( gamma(i) - z )';
    P_xz(:,:,i) = Wc * ( X_sig(i) - Xsa(1) ) * ( gamma(i) - z )';
end

P_z(:,:,1) = Wc0 * (gamma(1) - z)*(gamma(1) - z)';
P_z = sum(P_z, 3);

P_xz(:,:,1) = Wc * ( X_sig(1) - Xsa(1) ) * ( gamma(1) - z )';
P_xz = sum(P_xz, 3);



%% Kalman gain:

K = P_xz / P_z;

%keyboard

%% State and covariance estimate update:

X = X_s + K * (alpha - z);

X = X(1:2);


P = P_s - K * P_z * K';

P = P(1:2, 1:2);







%pause(0.1)

clc;
X
K
P_z
P_xz


%keyboard

