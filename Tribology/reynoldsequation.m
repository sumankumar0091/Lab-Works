clear all
close all

T1=3*h(i,k)^2;
T2= (p(i+1,k) - p(i+1,k))/(2*x);
T3= (h(i+1,k) - h(i+1,k))/(2*x);
T4=T1*T2*T3;
T5= (h(i,k)^3)/(x^2);
T6 = (p(i+1,k) - 2*p(i,k)+ p(i+1,k))/(z^2);
T7 = (2*l/w)^2;
T8 = (h(i,k+1) - h(i,k-1))/(2*z);
T9= (p(i,k+1) - p(i,k-1))/(z);
T10 = (h(i,k)^3)/(z^2);
T11= p(i,k+1) - p(i,k-1);
T12= 3*(h(i+1,k) - h(i-1,k))/(x);
T13=T4 + T5*T6;
T14 = 2*T5;
T15 = T7*T8*T1*T9;
T16=T7*T10*T11;
T17=2*T7*T10;
T18=T13 +T15 T16 - T12;
T19 = T14 + T17;
