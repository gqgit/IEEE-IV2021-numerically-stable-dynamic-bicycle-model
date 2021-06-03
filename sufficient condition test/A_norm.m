%Numerically Stable Dynamic Bicycle Model for Discrete-time Control (IEEE IV2021 workshop paper)
%this code checks the sufficient condition of numerical stability of our method with varying discretization steplength h.
%Author: Qiang Ge
%Date: 2021.06.03
%email: gq17@mails.tsinghua.edu.cn
u=0:0.1:15;
h1=0.01;
size = 15;
wid = 1;

norm2 = u;
for i = 1:length(u)
norm2(i)=norm(construct_A_hat(u(i), h1));
end
plot(u,norm2,'linewidth',wid)
%grid on
hold on

h2=0.05;
for i = 1:length(u)
norm2(i)=norm(construct_A_hat(u(i), h2));
end
plot(u,norm2,'linewidth',wid)
%grid on
hold on

h3=0.1;
for i = 1:length(u)
norm2(i)=norm(construct_A_hat(u(i), h3));
end
plot(u,norm2,'linewidth',wid)
%grid on

axis([0,max(u),0,1]);
hl = legend(strcat('$T_s$ = ',num2str(h1),' s'),strcat('$T_s$ = ',num2str(h2),' s'), strcat('$T_s$ = ',num2str(h3),' s'),'Location','Southeast');
set(hl,'Box','off','Interpret','latex','FontSize',size);
xlabel('longitudinal velocity $u_k$ (m/s)','Interpret','latex','FontSize',size)
ylabel('$\|\hat{A_k}\|_2$','Interpret','latex','FontSize',size)
yticks([0 0.5 1])
yticklabels({'0','0.5','1'})
set(gca,'FontName','Times New Roman','FontSize',14)
%ylabel('2-norm of $\hat A_k$','Interpret','latex','FontSize',size)


