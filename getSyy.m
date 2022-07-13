function [H,Hq,Cv,Fv,FI_v,FD_v] = getSyy(om_v,Sxx,A11n,A55n,A15n,B11n,B55n,B15n,Tc,para)


[M,K,Ms,xc,xcb,xi_v,aw,N,Ndof,Phi_x,V,D_om,Li_v,Di_v] = getMK(para);

[H,Hq,Cv,Fv,FI_v,FD_v] = solvefrf(M,K,Ms,V,D_om,Phi_x,om_v,Ndof,N,aw,xi_v,xc,Di_v,Li_v,A11n,A55n,A15n,B11n,B55n,B15n,Sxx,Tc,para);
                 
  

function [H,Hq,Cv,Fv,FI_v,FD_v]=solvefrf(M,K,Ms,V,D_om,Phi_x,om_v,Ndof,N,aw,xi_v,xc,Di_v,Li_v,A11n,A55n,A15n,B11n,B55n,B15n,Sxx,Tc,para)

Cd = para.Cd;
Cm = para.Ca + 1;
g = para.g; 
rho_f = para.rho_f;


    om_v = om_v(:)'; % make sure om_v is row vector 
   
    Nom=numel(om_v); 
    dom=om_v(2)-om_v(1); % delta omega 
    
    
    
    Coef_iter=0.5; % iteration steps <---------------------------
    Niter=50;
    
    
    % define intial vectors 

    % use deep water assumption here as d is big
    k_w = om_v.^2/g;  % wavenumber (based on linear wave theory) 
    
    u_i_v =  om_v.*aw.*exp(-xi_v*k_w);  % fluid velocity, matrix format, N x Nom 
    u_i_dot_v = om_v.^2.*aw.*exp(-xi_v*k_w); % fluid acceleration,  N x Nom 
    u_i_v(xi_v<0,:)=0; 
    u_i_dot_v(xi_v<0,:)=0;
    
    
    q_dot=zeros(Ndof,Nom); % q_dot=[alpha_dot theta_dot]' for rigid body case
    
    y_dot=q_dot(1)-(xi_v-xc)*q_dot(2); % transform back to physical coor 
    ur=u_i_v-y_dot;  % relative velocity for ith strip 
        
    
    % linearization coefficient 
    Sxx=Sxx.';
    Sur=abs(ur).^2.*Sxx;  
    gamma_old=sqrt(8/pi)*sqrt(sum(Sur*dom,2)); % linearised coefficient, initialise (sqrt(8/pi)*std(ur)) 

    
    % iteration to solve the equation because C is dependent on gamma
    for iter=1:Niter
        Cv=zeros(Ndof+1,Nom);
        Fv=zeros(Ndof,Nom);
        FI_v=zeros(Ndof,Nom); % inertia part
        FD_v=zeros(Ndof,Nom); % drag part
        
        H=zeros(N,Nom); %response spectrum,each row for each dof
        Hq=zeros(Ndof,Nom);
        
        ur=u_i_v-y_dot;  % relative velocity for ith strip 
        Sur=abs(ur).^2.*Sxx;  
        gamma_new=sqrt(8/pi)*sqrt(sum(Sur*dom,2)); % linearised coefficient, initialise (sqrt(8/pi)*std(ur))           

    
        gamma_T=gamma_old + ...
                  Coef_iter*(gamma_new-gamma_old); % iteration here
         
           
            % loop over omega   
            jj=0;
            for om=om_v % harmonic frequenccy from exp(i*om*t)
                
                jj=jj+1; 
                gamma=gamma_T; 

                u_v=u_i_v(:,jj);  % fluid velocity, vector of all strips
                u_dot_v=u_i_dot_v(:,jj); % fluid acceleration, vector of all strips 



                % Damping Matrix (needs update iteratively
                BDi_v=1/2*rho_f*Di_v.*Li_v.*Cd.*gamma; 
                BDi_v(xi_v<0)=0; % only effective below fluid surface 

                C=[sum(BDi_v) -sum(BDi_v.*(xi_v-xc));-sum(BDi_v.*(xi_v-xc)) sum(BDi_v.*(xi_v-xc).^2)] ;

                   
                % force vector (needs update iteratively
                F_I=pi/4*(Di_v).^2*rho_f*Cm.*Li_v.*u_dot_v;
                F_I(xi_v<0)=0; % only effective below fluid surface 
                
                F_D=1/2*rho_f*Cd.*Di_v.*Li_v.*u_v.*gamma;
                F_D(xi_v<0)=0; % only effective below fluid surface 
                
                FMi_v=F_I+F_D; 
                FMi_v(xi_v<0)=0; % only effective below fluid surface 

                F1=sum(FMi_v);
                F2=-sum((xi_v-xc).*FMi_v);
                F=[F1;F2];
                
                F_I_1=sum(F_I);
                F_I_2=-sum((xi_v-xc).*F_I);
                
                F_D_1=sum(F_D);
                F_D_2=-sum((xi_v-xc).*F_D);
                                
 
                % calculate response in frequency domain
                Ma=[A11n(jj) A15n(jj);A15n(jj) A55n(jj)];  % 2 by 2 matrix 
                Cp=[B11n(jj) B15n(jj);B15n(jj) B55n(jj)];  % 2 by 2 matrix 
                Ma=Tc.'*Ma*Tc;
                Cp=Tc.'*Cp*Tc;
                D_Y=-om^2*(Ms+Ma)+1i*om*(C+Cp)+K; % dynamic stiffness matrix 
                
                q_v=D_Y\F;  % generlised response (for rigid body case, these are the alpha & theta as the generalised coordinates)

                q_dot_v=1i*om*q_v; % generalised coor 
                
                y_dot(:,jj)=q_dot_v(1)-(xi_v-xc)*q_dot_v(2); % transform back to physical coor 
               
                H(:,jj)=q_v(1)-(xi_v-xc)*q_v(2); % y response (horizontal)           
                Hq(:,jj)=q_v; % y response (horizontal)

                Cv(:,jj)=[C(1,1);C(1,2);C(2,2)];
                Fv(:,jj)=F; 
                FI_v(:,jj)=[F_I_1;F_I_2]; 
                FD_v(:,jj)=[F_D_1;F_D_2];
                
                
            end
        
           gamma_old=gamma_new; % save old gamma before next iteration cycle   
          
    end
    
