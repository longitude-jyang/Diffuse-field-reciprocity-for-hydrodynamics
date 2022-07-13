% cal_disproots find the roots of the dispersion relation ktanh(kd)=om^2/g
% using Newton's method iteratively 
% the real roots and imaginary roots are found separately 

function [k0,kp]=cal_disproots(d,om)
% global g

g=9.81;
% iteration parameters  
nmax=50; 
epss=1; 
n=0; 


% for real root k0
    x=om^2/g+1; % initial value
    while epss>=1e-9 && n<=nmax 
        x1=x-(x*tanh(x*d)-om^2/g)/(tanh(x*d)+x*d*sech(x*d)^2); 
        epss=abs(x1-x);        
        x=x1;
        n=n+1; 
    end
    k0=x; 

    
% for complex roots kp
    C=d*om^2/g; % constant 

%     xp_v=om^2/g/d:pi/d:5; % vector of inital values (resolution of half wavelength equal the structure length)
    xp_v=om^2/g*linspace(d/100,2*d,100);
    
    
    Np=numel(xp_v);

    error=zeros(Np,1);
    kp=zeros(Np,1);


    epss=1; 
    n=0; 

    for jj=1:Np
        xp=xp_v(jj)*d;
        while epss>=1e-9&& n<=nmax 
            y=xp-(xp*tan(xp)+C)/(tan(xp)+xp*sec(xp)^2); 
            epss=abs(y-xp); 
            xp=y;n=n+1; 
        end
        n=0;
        epss=1;

        error(jj)=xp*tan(xp)+C;

        kp(jj)=xp/d;

    end
    
    kp=kp(kp>0);