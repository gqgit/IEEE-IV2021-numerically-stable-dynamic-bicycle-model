function [A] = construct_A_hat(u,h)
%Author: Qiang Ge
%Date: 2021.06.03
%email: gq17@mails.tsinghua.edu.cn
m = 1412;
Iz = 1536.7;
a = 1.06;
b = 1.85;
k1 = -128916;
k2 = -85944;

A = [       m*u/(m*u-h*(k1+k2))                  (h*(a*k1-b*k2)-h*m*u^2)/(m*u-h*(k1+k2));
         h*(a*k1-b*k2)/(Iz*u-h*(a^2*k1+b^2*k2))            Iz*u/(Iz*u-h*(a^2*k1+b^2*k2));  ];

end