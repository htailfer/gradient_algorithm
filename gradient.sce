function y=J(x)
    y=10*(x(2)-x(1)^2)^2-(x(1)-1)^2
endfunction
function y=nablaJ(x)
    y=([-40*x(1)*(x(2)-x(1)^2)+2*(x(1)-1);20*(x(2)-x(1)^2)])
endfunction
function y=hessJ(x)
    y=[-40*x(2)+120*x(1)^2+2,-40*x(1);-40*x(1),20]
endfunction
x=-1:0.01:2;y=-1:0.01:2;
z=zeros(length(x),length(y))
for i=1:length(x)
    for j=1:length(y)
      z(i,j)=J([x(i);y(j)])
    end
end
clf()
subplot(221)
surf(x,y,z)
rho=0.01
max_iter=2000
epsilon=1E-16



function[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho)
    nb_iter=0;
    u=u0;
    utot=u0;
    resid=nablaJ(u);   
    while(nb_iter<max_iter)&(norm(resid)>epsilon)
        u=u-rho*resid;
        utot=[utot,u];
        resid=nablaJ(u); 
        nb_iter=nb_iter+1;
    end
endfunction

u0=[1;0]
epsilon=1E-4
max_iter=20000
rho=0.01

[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
subplot(222)
contour(x,y,z,[0:0.5:2,3:10]);
plot(utot(1,:),utot(2,:));
title('rho=0,01')


subplot(223)

u0=[1;0]
epsilon=1E-4
max_iter=20000
rho=0.02
[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
contour(x,y,z,[0:0.5:2,3:10]);
plot(utot(1,:),utot(2,:));
title('rho=0,02')


subplot(224)

u0=[1;0]
epsilon=1E-4
max_iter=20000
rho=0.03
[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
contour(x,y,z,[0:0.5:2,3:10]);
plot(utot(1,:),utot(2,:));
title('rho=0,03')



function[utot,nb_iter]=NPF(nablaJ,hessJ,u0,epsilon,max_iter)
    nb_iter=0;
    u=u0;
    utot=u0;
    resid=nablaJ(u);
    while(nb_iter<max_iter)&(norm(resid)>epsilon)
        u=u-(inv(hessJ(u))*nablaJ(u))
        utot=[utot,u]
        resid=nablaJ(u)
        nb_iter=nb_iter+1 
    end
endfunction

u0=[1.1;0.1];
epsilon=1E-16;
max_iter=200;
rho=0.01;

[utot,nb_iter]=NPF(nablaJ,hessJ,u0,epsilon,max_iter)
err=utot-ones(utot);
normerr=sqrt(err(1,:).^2+err(2,:).^2)
plot(1:length(normerr),normerr,'k')
rho=0.01
max_iter=2000
[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
err1=utot-ones(utot);
normerr1=sqrt(err1(1,:).^2+err1(2,:).^2)
plot(1:length(normerr1),normerr1,'r')
rho=0.02
[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
err2=utot-ones(utot);
normerr2=sqrt(err2(1,:).^2+err2(2,:).^2)
plot(1:length(normerr2),normerr2,'g')
rho=0.03
[utot,nb_iter]=GPF2(J,nablaJ,u0,epsilon,max_iter,rho);
err3=utot-ones(utot);
normerr3=sqrt(err3(1,:).^2+err3(2,:).^2)
plot(1:length(normerr3),normerr3,'b')

title('Vitesse de convergence')

[u,nb_iter]=NPF(nablaJ,hessJ,u0,epsilon,max_iter)

