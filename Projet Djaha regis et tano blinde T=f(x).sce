//declaration des paramètres.
lambda=10;h=120;A=0.02;q=200;B=5;omega=0.5;
L=0.1;//épaisseur du mur en m
rho=2200;//mase volumique en Kg.m^-3
Cp=840;//capacité thermique J.m^-1.K^-1
T0=30;Tf=30;T2=30;T3=30;//température en degré celcius
t=180;//est pris de façon arbitraire t=60 ,t=120,t=180,t=360 et t=600t=180;//est pris de façon arbitraire t=60 ,t=120,t=180,t=360 et t=600
tfinal=180;dt=0.5;m=tfinal/t//temps en seconde


//calcul des differents coefficients des matrices.
p1=6*dt/(rho*Cp*A*L);p2=dt*(q*A*L/2-(rho*Cp*A*B*L*omega*cos(omega*t))/6+(lambda*A/L)*(T0+B*sin(omega*t)));p3=A*h*Tf*dt;

//construction des matrices
MatA=[4 1;1 2]
InvMatA=inv(MatA);

MatB(1,1)=2*dt*lambda*A/L;
MatB(1,2)=-dt*lambda*A/L;
MatB(2,1)=-dt*lambda*A/L;
MatB(2,2)=(dt*lambda*A+A*L*h)/L;

VectF(1)=p2;VectF(2)=p3;

//resolution du systeme d'equations
t=0;Tn(1)=T2;Tn(2)=T3;
for n=1:m
    t=t+dt;
    Tnplusun=Tn+p1*InvMatA*(-MatB*Tn+VectF);
    Tnsolution(n,1)=Tnplusun(1,1);
    Tnsolution(n,2)=Tnplusun(2,1);

T2=[T2;Tnsolution(:,1)];
T3=[T3;Tnsolution(:,2)];
t=0:dt:tfinal;
disp(t,Tnplusun)//afficher toutes les valeurs de Tnplusun 
//Pour chaque figure fixons T1=T0+B*sin(omega*t) en faisant varié t
//traçons chaque courbe en fonction de l'espace en fixant x=0:0.1:0.2 et et trouvons chaque valeur de la temperature pour t fixé par un encardrement décroissante
end
