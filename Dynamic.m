function [x1_dot,x2_dot] = Dynamics(x1,x2,u)
%#codegen

theta1 = 0.1;
theta2 = 0.1;
theta3 = -0.2;



x1_dot = theta1*x1^2 + x2;

x2_dot = theta2*x1*x2 + theta3*x1 + (1 + x1^2)*u;