function [e1,e2,u,ac_dot] = Controller(t,x1,x2,ac)
%#codegen

theta1 = 0.1;
theta2 = 0.1;
theta3 = -0.2;

kc1 = 1.0;
kc2 = 2.5;

k1 = 5;
k2 = 5;

x1_dot = theta1*x1^2 + x2;
%x2_dot = theta2*x1*x2 + theta3*x1 + (1 + x1^2)*u;

yd = 0.2*sin(2.5*t);
yd_dot = 0.2*2.5*cos(2.5*t);
yd_ddot = -0.2*2.5^2*sin(2.5*t);

%WORKING
% yd = 0.2;
% yd_dot = 0;
% yd_ddot = 0;

z1 = x1-yd;
z1_dot = x1_dot-yd_dot;

%tau1 = (kc1/(2*z1))*log((kc1+z1+yd)*(kc1-yd)/((kc1-z1-yd)*(kc1+yd)));
if abs(z1) == 0
    tau1 = kc1^2/(kc1^2-yd^2);
    tau1_dot_dz = kc1/2*(1/(kc1+z1+yd)^2 + 1/(kc1-z1-yd)^2);
else    
    tau1 = (kc1/(2*z1))*log((kc1+z1+yd)*(kc1-yd)/((kc1-z1-yd)*(kc1+yd)));
    tau1_dot_dz = (1/z1)*(kc1^2/(kc1^2-(z1+yd)^2)-tau1)*z1_dot;
end


tau1_dot = tau1_dot_dz + yd_dot*kc1^2*(z1+2*yd)/(kc1^2-(z1+yd)^2)/(kc1^2-yd^2);

a1 = -theta1*(x1^2)-k1*z1+yd_dot*(kc1^2-x1^2)*tau1/(kc1^2);
z2 = x2 - a1;

%tau2 = (kc2/(2*z2))*log((kc2+z2+a1)*(kc2-a1)/((kc2-z2-a1)*(kc2+a1)));
if abs(z2) == 0
   tau2 = kc2^2/(kc2^2-a1^2);
else
    tau2 = (kc2/(2*z2))*log((kc2+z2+a1)*(kc2-a1)/((kc2-z2-a1)*(kc2+a1)));
end

%a1_dot is not used when Dynamic Surface Approximation done
a1_dot = -theta1*(2*x1*x1_dot)-k1*z1_dot + yd_dot*(kc1^2-x1^2)*tau1_dot/(kc1^2) - 2*yd_dot*x1*x1_dot*tau1/(kc1^2) + yd_ddot*(kc1^2-x1^2)*tau1/(kc1^2);
%comment ac_dot and replace ac with a1 for NO Dynamic Surface Approx.
gamma = 0.02;
ac_dot = 1/gamma*(-ac + a1);
u = 1/(1+x1^2)*(-theta3*x1-theta2*x1*x2-k2*z2+ac_dot*tau2*(kc2^2-x2^2)/kc2^2-(kc2^2-x2^2)*kc1^2*z1/(kc1^2-x1^2)/kc2^2);


e1 = x1 - yd;
e2 = x2 - yd_dot;