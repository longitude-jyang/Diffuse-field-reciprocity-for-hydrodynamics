function [M,K,Ms,xc,xcb,xi_v,aw,N,Ndof,Phi_x,V,D_om,Li_v,Di_v]=getMK(para)


    rho_s = para.rho_s; 
    rho_f = para.rho_f;
    s_alpha = para.s_alpha;
    Ca = para.Ca;
    D1 = para.D1;
    D2 = para.D2;
    L1 = para.L1;
    L2 = para.L2;
    g  = para.g; 

    aw = 1;
    %-------------------------------------------------------------------------
        %  Design Parameter for the cylinder 

        %  total length of the cylinder 
            L = L1+L2; 

        %  Section -  Diameters
            deltaD =0.02*5;
            Dae = D1;
            Dbe = D2;
            Dce = D2;
        %  Section -  Density    
            rho_a=rho_s;
            rho_b=rho_s;
            rho_c=rho_s*45;%

        %  Section - Lengths    
            La=L1+6;
            Lb=L2-6-10;
            Lc=10;


            Na=La;
            Nb=Lb;
            Nc=Lc;
            xi_v=[linspace(-L1,La-L1,Na).';...
                linspace(La-L1+Lb/Nb,L2-Lc,Nb).';...
                linspace(L2-Lc+Lc/Nc,L2,Nc).'] ;


            N=Na+Nb+Nc ;% number of strips 

        %-------------------------------------------------------------------------
        %  parameters for strips 


        % loop over the strips 
        Mi_v=zeros(N,1); % vector for mass of all strips 
        Mia_v=zeros(N,1); % vector for added mass (the ones above fluid surface will be zeros
        Ii_v=zeros(N,1); % vector for 2nd moment of area
        Li_v=zeros(N,1) ;

        Di_v=zeros(N,1) ;
        Di_v_I=zeros(N,1) ;


        for ii=1:N 

                if ii<=Na
                    rho_i=rho_a;
                    Di=Dae;
                    Di_I=Dae-deltaD;
                    Li=La/Na;
                elseif ii>Na && ii<=Na+Nb
                    rho_i=rho_b;
                    Di=Dbe;
                    Di_I=Dbe-deltaD;
                    Li=Lb/Nb;
                elseif ii>Na+Nb 
                    rho_i=rho_c;
                    Di=Dce;
                    Di_I=Dce-deltaD;
                    Li=Lc/Nc;
                end
                Di_v(ii)=Di; % diameter of ith strip, outer
                Di_v_I(ii)=Di_I; % inner diameter
                Vi=pi*(Di/2)^2*Li; % volume of ith strip 
                Vi_I=pi*(Di_I/2)^2*Li;

                Mi=rho_i*(Vi-Vi_I) ; % mass of ith strip

                Ii=pi*((Di/2)^4-(Di_I/2)^4)/4; % 2nd moment of area, for circular cross-section pi*radius^4/4

                Bi=Vi*rho_f; % buoyancy of ith trip ( equal to the mass of displaced fluid)  
                Mia=Ca*Bi; % added mass from fluid, Bi is buoyancy of ith trip 

                % put them into vector form 
                Mi_v(ii)=Mi;
                Mia_v(ii)=Mia; 

                Ii_v(ii)=Ii;
                Li_v(ii)=Li;

        end 

        Mia_v(xi_v<=0)=0; % no fluid above fluid surface 

        xc = sum(Mi_v.*xi_v)/sum(Mi_v); % get centre of gravity 
        xm = xi_v(end-Nc); % coor for spring 

        xcb=sum(Mia_v.*xi_v)/sum(Mia_v); % get centre of buoyancy 

        Mt = sum(Mi_v); % total structure mass
        Mat = sum(Mia_v);
        
        
        % -----------------------------------------------
                Ndof=2; 
% 
                Phi_x=[ones(N,1)/L (xi_v-xc)/L];% size N x 2
     
                M=zeros(Ndof,Ndof); 
                Ms=zeros(Ndof,Ndof);
                K=zeros(Ndof,Ndof); 

                % mass matrix
                M(1,1)=sum(Mi_v+Mia_v); 
                M(2,2)=sum((Mi_v+Mia_v).*(xi_v-xc).^2);
                
                Ms(1,1)=sum(Mi_v); 
                Ms(2,2)=sum((Mi_v).*(xi_v-xc).^2);

                % stiffness matrix
                K(1,1)=s_alpha;
                K(2,2)=sum((Mi_v-Mia_v/Ca)*g.*(xi_v-xc))+s_alpha*((xm-xc).^2);
                K(1,2)=-s_alpha*(xm-xc);
                K(2,1)=K(1,2); 
%                 
                % undamped natural frequency 
                
                [V,D_om]=eig(K,M); 
                D_om=diag(D_om);
                [sorter,index]=sort(D_om);
                V=V(:,index);
                D_om=D_om(index);
                
                om_n=sqrt(D_om);
                
                % Mass normalization 
                for ii=1:Ndof        
                    V(:,ii)=V(:,ii)/sqrt((V(:,ii).'*M*V(:,ii)));
                end
                

       