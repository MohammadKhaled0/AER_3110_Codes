% Mohammed Khalid Gamal Ali / sec:2 / B.N:14 
% Solution of the flow past Joukowski airfoil
close all; clear all ; clc

%% Input data
V_inf = 100;                                    % Free stream velocity ( m/s)
Cord = 1;                                       % Chord length of the airfoil (m)
Max_Camber = 7 / 100;                           % max. camber to chord ratio
Max_Thiccness = 6 / 100;                        % max. thickness to chord ratio
A_O_A = 5;                                      % Angle of attack (deg)
Max_No_Of_Upper_points = 101;                   % Total number of points on the airfoil upper surface
Max_No_Of_Lower_points = 101;                   % Total number of points on the airfoil lower surface
%% Calculation of constants appearing in the equations
Sin_a = sind(A_O_A); 
Cos_a = cosd(A_O_A);
b = Cord / 4; e = Max_Thiccness / 1.3;          %Transformation Constant
Beta = 2 * Max_Camber; 
Cos_b = cos(Beta); 
Sin_b = sin(Beta);
a = b * (1 + e) / Cos_b; 
X_O = -b * e;
Y_O = a * Sin_b;
%% Calculation of the x & y coordinates and the velocity on the airfoil upper surface
for i=1:1:Max_No_Of_Upper_points
Theta = ( pi) * (i - 1) / (Max_No_Of_Upper_points - 1);
if (i == 1) ; Theta = .01 * pi / 180; end
Cos_theta = cos(Theta); 
Sin_theta = sin(Theta); 
X_1 = 2 * b * Cos_theta;
Y_1 = 2 * b * e * (1 - Cos_theta) * Sin_theta + 2 * b * Beta * Sin_theta^2;
X_Upper(i)=X_1; 
Y_Upper(i)=Y_1;
r = b * (1 + e * (1 - Cos_theta) + Beta * Sin_theta); 
X= r * Cos_theta; 
Y = r * Sin_theta;
X_dash = X - X_O; Y_dash = Y - Y_O; rd = sqrt(X_dash * X_dash + Y_dash * Y_dash);
cos_theta_dash = X_dash / rd; sin_theta_dash_A = Y_dash / rd;
cos_theta_dash_A = cos_theta_dash * Cos_a + sin_theta_dash_A * Sin_a; 
sin_theta_dash_A = sin_theta_dash_A * Cos_a - cos_theta_dash * Sin_a;
sinapb = Sin_a * Cos_b + Cos_a * Sin_b;
vrd = V_inf * (1 - (a * a / (rd * rd))) * cos_theta_dash_A;
vtd = -V_inf * ((1 + (a * a / (rd * rd))) * sin_theta_dash_A + 2 * (a / rd) * sinapb);
cos2t = Cos_theta * Cos_theta - Sin_theta * Sin_theta; sin2t = Sin_theta * Cos_theta + Cos_theta * Sin_theta;
TA = vrd * cos_theta_dash - vtd * sin_theta_dash_A; TB = -(vrd * sin_theta_dash_A + vtd * cos_theta_dash);
TC = 1 - ((b / r) ^ 2) * cos2t; TD = ((b / r) ^ 2) * sin2t;
V2 = (TA * TA + TB * TB) / (TC * TC + TD * TD); V = sqrt(V2); Va_Upper(i) =V/V_inf;
end
for i=1:1:Max_No_Of_Lower_points
Theta = pi + pi* (i - 1) / (Max_No_Of_Lower_points - 1);
if (i == Max_No_Of_Lower_points) Theta = 2 * pi - .01 * pi / 180; end
Cos_theta = cos(Theta); 
Sin_theta = sin(Theta); 
X_1 = 2 * b * Cos_theta;
Y_1 = 2 * b * e * (1 - Cos_theta) * Sin_theta + 2 * b * Beta * Sin_theta^2;
X_Lower(i)=X_1; 
Y_Lower(i)=Y_1;
r = b * (1 + e * (1 - Cos_theta) + Beta * Sin_theta); 
X= r * Cos_theta; 
Y = r * Sin_theta;
X_dash = X - X_O; Y_dash = Y - Y_O; rd = sqrt(X_dash * X_dash + Y_dash * Y_dash);
cos_theta_dash = X_dash / rd; sin_theta_dash_A = Y_dash / rd;
cos_theta_dash_A = cos_theta_dash * Cos_a + sin_theta_dash_A * Sin_a; 
sin_theta_dash_A = sin_theta_dash_A * Cos_a - cos_theta_dash * Sin_a;
sinapb = Sin_a * Cos_b + Cos_a * Sin_b;
vrd = V_inf * (1 - (a * a / (rd * rd))) * cos_theta_dash_A;
vtd = -V_inf * ((1 + (a * a / (rd * rd))) * sin_theta_dash_A + 2 * (a / rd) * sinapb);
cos2t = Cos_theta * Cos_theta - Sin_theta * Sin_theta; sin2t = Sin_theta * Cos_theta + Cos_theta * Sin_theta;
TA = vrd * cos_theta_dash - vtd * sin_theta_dash_A; TB = -(vrd * sin_theta_dash_A + vtd * cos_theta_dash);
TC = 1 - ((b / r) ^ 2) * cos2t; TD = ((b / r) ^ 2) * sin2t;
V2 = (TA * TA + TB * TB) / (TC * TC + TD * TD); V = sqrt(V2); Va_Lower(i) =V/V_inf;
end
figure
subplot(4,1,4);
plot (X_Lower,Y_Lower); hold on; grid on
plot (X_Upper,Y_Upper);
axis equal;
xlabel('airfoil X-axis')
ylabel('airfoil Y-axis')
title(['Joukwoski Airfoil "camb.ratio"=',num2str((Max_Camber)),',"thicc.ratio"=',num2str((Max_Thiccness)),',"angle of attack"=',num2str((A_O_A)),''])
subplot(4,1,1:3);
plot (X_Lower,Va_Lower); hold on; grid on
plot (X_Upper,Va_Upper);
xlabel('airfoil X-axis')
ylabel('Velocity ratio "V/Vinf"')
title(['3- Velocity distribution over the airfoil surface for angle of attack 5^{\circ}'])
