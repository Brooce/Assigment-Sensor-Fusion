

if (k < 4)
    return
end


%% "iso range ratio" circle centre and radius calculation based on alpha

%
[Ck, Rk] = get_geo_data(x_vec_all(1,:), x_vec_all(k,:), P_r_filt_ratio(k,1));
[Ck_1, Rk_1] = get_geo_data(x_vec_all(1,:), x_vec_all(k-1,:), P_r_filt_ratio(k-1,1));
[Ck_2, Rk_2] = get_geo_data(x_vec_all(1,:), x_vec_all(k-2,:), P_r_filt_ratio(k-2,1));
%}

%
x1 = Ck(1);
y1 = Ck(2);
r1 = Rk;
x2 = Ck_1(1);
y2 = Ck_1(2);
r2 = Rk_1;
%}

%{
x1 = 1;
y1 = 1;
r1 = 1;
x2 = 2;
y2 = 2;
r2 = 2;
%}

%% Symbolic equation to find the crossing points

syms x y

circle1 = [(x - x1)^2 + (y - y1)^2 == r1];
circle2 = [(x - x2)^2 + (y - y2)^2 == r2];


crossing_points = solve(circle1, circle2, 'x','y');

%cross_x = double(odp.x);
%cross_y = double(odp.y);

cross1 = [ double(crossing_points.x(1))     double(crossing_points.y(1))  ];
cross2 = [ double(crossing_points.x(2))     double(crossing_points.y(2))  ];

%crossing_points = [ double(odp.x)'    ; double(odp.y)'  ];


%% Choosing one of these points basing on the third circle

distsq_Ck2_c1 = (cross1 - Ck_2) * (cross1 - Ck_2)';
distsq_Ck2_c2 = (cross2 - Ck_2) * (cross2 - Ck_2)';

if ( distsq_Ck2_c1 - Rk_2 < distsq_Ck2_c2 - Rk_2 )
    x_state_prev = cross1';
else
    x_state_prev = cross2';
end



