function [c,ceq] = constraintA(x)
global minCL_Const;
global maxCM_Const
global LD_d;
global CL_0_d;
global CM_Max_d;

x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
x5=x(5);


model_LD= [1,x1, x2, x3, x4, x5, x1.^2, x2.^2, x3.^2, x4.^2, x5.^2,x1.*x2,x2*x3,x3*x4,x4*x5]*LD_d;
model_CL0=[1,x1, x2, x3, x4, x5, x1.^2, x2.^2, x3.^2, x4.^2, x5.^2,x1.*x2,x2*x3,x3*x4,x4*x5]*CL_0_d;
model_CM_Max=[1,x1, x2, x3, x4, x5, x1.^2, x2.^2, x3.^2, x4.^2, x5.^2,x1.*x2,x2*x3,x3*x4,x4*x5]*CM_Max_d;




c=[minCL_Const-model_CL0;model_CM_Max-maxCM_Const];
ceq=[];
end

