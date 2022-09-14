% Mohammed Khalid Gamal Ali / sec:2 / B.N:14
% Solution of the flow past Joukowski airfoil
close all; clear all ; clc ; digits(4);
%% Input data
V_inf = 100;                % Free stream velocity ( m/s) % Airfoil Geometric Characteristics
cord = 1 ;                  % Chord length of the airfoil (m)
Max_Thiccness = 6/100;      % max. thickness to chord ratio
Max_Camber = 7/100;         % max. camber to chord ratio
%% Calculation of constants appearing in the equations
Theta = linspace(1,360-1,2*180)*pi/180;
b = cord/4; e = Max_Thiccness/1.3;
beta = 2 * Max_Camber;
a = b * (1+e);
x_0 = -b * e;
y_0 = a * beta;
%% Calculation of the x & y coordinates and the velocity on the airfoil
X = -2.5*(b+x_0) :b/100:2.5*(b+x_0);
Y = -1.25*(a+y_0):a/100:1.25*(a+y_0);
x_1 = 2*b*cos(Theta);
y_1 = 2*b*e*(1-cos(Theta)).*sin(Theta)+2*b*beta*sin(Theta).^2;
r = b*(1+e*(1-cos(Theta))+beta*sin(Theta));
Z_1 = x_1+1i*y_1;
x = r.*cos(Theta);
y = r.*sin(Theta);
Z_D = x + y*1i;
z0 = x_0 + y_0*1i; Z_D = Z_D-z0;
r_D = abs(Z_D);
Theta_D= angle(Z_D);
A_O_A = 5 *pi/180;
V_r_D = V_inf* (1-(a^2./r_D.^2)).*cos(Theta_D-A_O_A); V_teta_dash = -V_inf*((1+(a^2./r_D.^2)).* sin(Theta_D-A_O_A)+2*(a./r_D)*sin(A_O_A+beta));
A = V_r_D.*cos(Theta_D)- V_teta_dash.* sin(Theta_D);
B = -V_r_D.*sin(Theta_D)- V_teta_dash.* cos(Theta_D);
C = 1 - b^2 ./ r.^2 .* cos(2*Theta); D = b^2 ./ r.^2 .* sin(2*Theta);
dW_dZ1=(A+B.*1i)./(C+D.*1i); u = real(dW_dZ1); v = imag(dW_dZ1); V1 = abs(dW_dZ1);
C_p = 1-(V1/V_inf).^2;
%% A Drawing of the airfoil with dimensions. & Velocity and Pressure coeff. distribution over the airfoil surface for angle of attack 5
figure(2)
subplot(4,1,1:3); hold on ;grid on;axis tight
plot(x_1,V1/V_inf,'color',[25/255, 150/255, 145/255],'LineWidth',1.75)
figure(3)
subplot(4,1,1:3); hold on ;grid on;axis tight
plot(x_1,C_p,'color',[25/255, 150/255, 145/255],'LineWidth',1.75)
Cl = 1 / cord * trapz( - C_p , x_1 );

A_O_A=linspace(-4,13,360)*pi/180;
for ii=1:length(A_O_A)
    V_r_D = V_inf* (1-(a^2./r_D.^2)).*cos(Theta_D-A_O_A(ii)); V_teta_dash = -V_inf*((1+(a^2./r_D.^2)).* sin(Theta_D-A_O_A(ii))+2*(a./r_D)*sin(A_O_A(ii)+beta));
    A = V_r_D.*cos(Theta_D)- V_teta_dash.* sin(Theta_D); B = -V_r_D.*sin(Theta_D)- V_teta_dash.* cos(Theta_D); C = 1 - b^2 ./ r.^2 .* cos(2*Theta); D = b^2 ./ r.^2 .* sin(2*Theta);
    dW_dZ1= (A+B.*1i)./(C+D.*1i); V1 = abs(dW_dZ1);
    C_p = 1-(V1/V_inf).^2;
    Cl(ii) = 1 / cord * trapz( - C_p , x_1 );
    CD(ii) = trapz( - C_p , y_1);
    C_m_LE(ii)=1/cord*(trapz(-C_p.*(x_1+cord/2),x_1+cord/2) + trapz(-C_p.*y_1,y_1) );
    X_cp(ii) = trapz(-C_p.*x_1,x_1)./Cl(ii)/cord ;
end
alpha = A_O_A;
C_l = 2*pi*(1+e)*sin(beta+alpha);
%% Joukwoski Airfoil
figure(1)
title(['Joukwoski Airfoil camb.ratio = ',num2str((Max_Camber)),',thicc.ratio = ',num2str((Max_Thiccness)),', \alpha = 5',' \circ'])
hold on ;grid on;axis equal;
ylabel('Y')
xlabel('X')
P1 = plot(Z_1,'LineWidth', 1.75,'color',[25/255, 150/255, 145/255]);

%% streamlines around Jouckowsky airfoil
r_D=linspace(a,3*a,360);
Theta_D=linspace(0,2*pi,360); alpha=pi/180*5; Gamma=4*pi*V_inf*a*sin(alpha+beta);
[Theta_D,R_dash]= meshgrid(Theta_D,r_D); Z_D=R_dash.*exp(1i*Theta_D)+x_0+1i*y_0;
Z_1=Z_D+b^2./Z_D;
W =  V_inf   .*exp(-1i.*alpha(length(alpha))).*Z_D...
    + V_inf*a^2*exp( 1i *alpha(length(alpha)))./Z_D...
    + 1i*log(Z_D)*Gamma/2/pi; psi=imag(W);
psi1=psi+b^2./psi;

figure
title('streamlines around Jouckowsky airfoil')
hold on
axis equal
plot(Z_1(1,:),'color',[51/255, 153/255, 255/255]);
fill(real(Z_1(1,:)),imag(Z_1(1,:)),[51/255, 153/255, 255/255])
box on;
contour(real(Z_1),imag(Z_1),psi1,'levelstep',3,'color',[95/255, 4/255, 150/255]); xlim([-2/3*cord,2/3*cord]);
ylim([-1.5*a,1.5*a]);
%% 
figure(2)
subplot(4,1,1:3); title('V/V_\infty along airfoil')
ylabel('V/V_\infty')
legend('\alpha =5')
subplot(4,1,4); hold on ;grid on; axis equal
plot(x_1,y_1,'color',[25/255, 150/255, 145/255],'LineWidth', 1.75)
xlabel('X')
%%
figure(3)
subplot(4,1,1:3); title('\fontsize{10}C_p along airfoil')
ylabel('\fontsize{15}C_p')
legend('\alpha =5')
legend('Location','east')
subplot(4,1,4); hold on ;grid on; axis equal
plot(x_1,y_1,'color',[25/255, 150/255, 145/255],'LineWidth', 1.75)
xlabel('airfoil X-axis')
%% The variation of the lift coefficient with the angle of attack.
AoA=linspace(-5,10,length(Cl))*pi/180;
figure(4)
subplot(2,1,1)
plot(AoA*180/pi,Cl,'color',[25/255, 150/255, 145/255],'LineWidth', 1.75)
title('C_L Vs \alpha')
hold on
grid on;axis tight
xlabel('\alpha')
ylabel('C_L')
%% The variation of the moment coefficient with the angle of attack
subplot(2,1,2)
plot(AoA*180/pi,C_m_LE,'color',[25/255, 150/255, 145/255],'LineWidth', 2)
title('C_M_L_E Vs \alpha')
xlabel('\alpha')
ylabel('C_m')
grid on;axis tight


