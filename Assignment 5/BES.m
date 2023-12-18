% Abhishek Ghosh ME21BTECH11001
% ME3030 Assignment 5

% function for zdot
function zdot=BES(t,z)
global m1 m2 J1 J2 a b g A omega

%initial conditions & derivatives
xc1=1;
yc1=1;
xc2=1;
yc2=1;
xcd1=0;
ycd1=0;
xcdd1=0;
ycdd1=0;
xcd2=0;
ycd2=0;
xcdd2=0;
ycdd2=0;

% Mass matrix
M=diag([m1 m1 J1 m2 m2 J2]);
% Force matrix
F=[0 -m1*g 0 0 -m2*g 0]';
x1=z(1);
y1=z(2);
theta1=z(3);

x2=z(4);
y2=z(5);
theta2=z(6);

x1d=z(7);
y1d=z(8);
theta1d=z(9);

x2d=z(10);
y2d=z(11);
theta2d=z(12);
% Uqdd=v

% u-> coeeficient of qdd matrix

U=[1 0 a*sin(theta1)+b*cos(theta1) -1 0 a*sin(theta2)+b*cos(theta2);
   0 1 b*sin(theta1)-a*cos(theta1) 0 -1 b*sin(theta2)-a*cos(theta2);
   1 0 -a*sin(theta1)-b*cos(theta1) 0 0 0;
   0 1 a*cos(theta1)-b*sin(theta1) 0 0 0];

% v-> independent of qdd terms

v=[theta1d^2*(b*sin(theta1)-a*cos(theta1)) + theta2d^2*(b*sin(theta2)-a*cos(theta2));
   theta1d^2*(-a*sin(theta1)-b*cos(theta1)) + theta2d^2*(-a*sin(theta2)-b*cos(theta2));
   theta1d^2*(a*cos(theta1)-b*sin(theta1));
   theta1d^2*(b*cos(theta1)+a*(sin(theta1)))];

% acc-> derivatives of x1 y1 theta1 x2 y2 theta2

acc=M\F+(M^(-0.5))*pinv(U*(M^(-0.5)))*(v-U*(M\F));

zdot=[z(7) z(8) z(9) z(10) z(11) z(12) acc']';