% Mohammad Khaled Gamal Ali / sec:2 / B.N:50 
% LAB Assignment / Submitted to: Eng.Karim Ahmed
clear all
clc
%% Input data x
x_Lower=[1.00000,  0.99000,  0.98000,  0.96672,  0.94908,  0.92565,  0.89453,  0.85321,  0.79832,  0.72542,  0.62860,  0.50000,  0.37855,  0.28653,  0.21680,  0.16397,  0.12394,  0.09360,  0.07061,  0.05320,  0.04000,  0.03000,  0.02000,  0.01000,  0.00000];
x_Upper=[0.00000,  0.01000,  0.02000,  0.03000,  0.04000,  0.05320,  0.07061,  0.09360,  0.12394,  0.16397,  0.21680,  0.28653,  0.37855,  0.50000,  0.62860,  0.72542,  0.79832,  0.85321,  0.89453,  0.92565,  0.94908,  0.96672,  0.98000,  0.99000,  1.00000];
x = [x_Lower, x_Upper];

%% Y from the x data
y_Upper=+0.594689181*(0.298222773*sqrt(x_Upper)-0.127125232*x_Upper-0.357907906*x_Upper.^2+0.291984971*x_Upper.^3-0.105174606*x_Upper.^4);
y_Lower=-0.594689181*(0.298222773*sqrt(x_Lower)-0.127125232*x_Lower-0.357907906*x_Lower.^2+0.291984971*x_Lower.^3-0.105174606*x_Lower.^4);
y = [y_Lower, y_Upper];

%% Input data Cp
Cp_Upper= [-0.83000,  0.53900,  0.80000,  0.99470,  1.05000,  1.07000,  1.05560,  1.04170,  0.98080,  0.89740,  0.79500,  0.67130,  0.54330,  0.38950,  0.24340,  0.15300,  0.08340,  0.01400, -0.04860, -0.09040, -0.15990, -0.20170, -0.26430, -0.35600, -0.4000];
Cp_Lower= [-0.40000, -0.35510, -0.27826, -0.22260, -0.18087, -0.13913, -0.09739, -0.04350,  0.00695,  0.04174,  0.11130,  0.18080,  0.23650,  0.27100,  0.28500,  0.26440,  0.23650,  0.17390,  0.09040,  0.00000, -0.12520, -0.18000, -0.31650, -0.51130, -0.83000];
Cp=[Cp_Lower, Cp_Upper];
%% A.O.A.
a=2;
%% I.C. at point x(i+0.5) where i+0.5 is the midpoint of i and i+1
 
C=1;
C_Normal=0;
C_Axial=0;
Cm=0;
%% Calculations of  normal and axial forces
for n = 1:length(x)-1 
    
delta_xi = x(n+1)-x(n);
delta_yi = y(n+1)-y(n);

xi_Mid_Point = 0.5*(x(n)+x(n+1));
yi_Mid_Point = 0.5*(y(n)+y(n+1));
Cpi_Mid_Point = 0.5*(Cp(n)+Cp(n+1));

C_Normal =C_Normal  + Cpi_Mid_Point * delta_xi;
C_Axial =C_Axial  -Cpi_Mid_Point * delta_yi;

Cm =Cm + ( -(Cpi_Mid_Point * delta_xi)*xi_Mid_Point -(Cpi_Mid_Point * delta_yi)*yi_Mid_Point ); 

end
%% getting CL,CD and Cm from the normal and axial forces
CL = C_Normal*cos(a*pi/180) - C_Axial*sin(a*pi/180);
CD = C_Axial*cos(a*pi/180) - C_Normal*sin(a*pi/180);
%% Printing the Values of CL,CD,Cm
disp('CL=')
disp(CL)
disp('CD=')
disp(CD)
disp('Cm=')
disp(Cm)
%% Plotting of the X/c vs Cp
plot(x/C,Cp,'p','Color',[0.5,0.2,0.9])
hold on 
plot(x/C,Cp,'g')
xlabel('x/c')
ylabel('Cp')
legend( "This Solution","Stram Fun. Solution" )
str = {'CL=',CL,'CD=',CD,'Cm=',Cm};

text(0.2,-0.5,str,'FontSize',8)