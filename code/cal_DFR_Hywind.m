% cal_DFR_Hywind.m  uses Hywind example to demontrate the DFR application,
% once the hydrodynamic coefficients are known, where a spreading sea is considered
% and compares it with calucation based on Morison's equation, 
% where a uni-directional wave is assumed

% 17/12/2019 @ JD1, Cambridge  [J Yang] 

% --------------------------------------------
% (0)  parameters 
% --------------------------------------------
    
   % general 
   para.g = 9.81;

   % Morison's
   para.Ca = 1; 
   para.Cd = 1; 

   % water 
   para.d = 300;  % water depth
   para.rho_f = 1.025e3; % density 

   % structure 
   para.rho_s = 8.5e3; % density  (steel )
   para.L1 = 87;  % length above water 
   para.L2 = 120; % length below water 
 
   para.D1 = 6.5; % diameter (narrow section)
   para.D2 = 9.4; % diameter (fat section) 

   % mooring 
   para.s_alpha = 3.8e9;  % spring stiffness in the alpha direction

   % define a common frequency vector
   Nom = 50;
   om_v = linspace(0.4,1.4,Nom);
   dom  = om_v(2)-om_v(1);

   % random wave spectrum 
   [Sxx, ~, ~ ] = jonswap(om_v, 2.5, 2*pi/0.8); % jonswap wave 

% ------------------------------------------------------
%  (1) read presaved added mass and damping coefficients
% --------------------------------------------
    load('Hywind_AddMassDamping.mat');

    A11n = interp1(A11(:,1),A11(:,2),om_v); % A11 = A22
    A55n = interp1(A44(:,1),A44(:,2),om_v); % A44 = A55
    A15n = -interp1(A24(:,1),A24(:,2),om_v); % A24 = A15

    B11n = interp1(B11(:,1),B11(:,2),om_v);%

    [B44s,ia] = unique(B44(:,1));
    B55n = interp1(B44s,B44(ia,2),om_v);%
    B15n = -interp1(B24(:,1),B24(:,2),om_v);%   


% --------------------------------------------
% (2)  construct dynamic model of Hywind 
% --------------------------------------------   
 
    % M and K matrix from FWT code via Lagrange forumation with strips 
    [M,K,Ms,xc,xcb] = getMK(para);  
   
    % transformation matrix to tranform coordinate    
    % the added mass/damping was with respsect to SWL, 
    % but the response of FWT we calculate is with respect to center of mass   

    Lc = xc; % centre of mass below SWL  (if Lc=xcb, the added mass matrix should be diagonalised)
    Tc = [1 Lc;0 1]; % coordinate transformation matrix    
  
    % get FWT response transfer function from Morison's equation for a unit
    % wave amplitdue aw=1 
    [H,Hq,Cv,Fv,FI_v,FD_v] = getSyy(om_v,Sxx,A11n,A55n,A15n,B11n,B55n,B15n,Tc,para);


% --------------------------------------------
%  (3)  caculate blocked force and response from DFR and Morison's and compare 
% --------------------------------------------    
    % intialise vectors 
    Sff_om=zeros(Nom,3); % force spectrum
    Szz_om=zeros(Nom,3); % response cross spectrum
    
    SyyQ_om=zeros(Nom,3); % response cross spectrum from Morison's equation 
    SyyQ2_om=zeros(Nom,3); % response cross spectrum from Morison's equation 
    
    SffQ_om=zeros(Nom,3);
    SffQ_I_om=zeros(Nom,3);
    SffQ_D_om=zeros(Nom,3);
    
    E2n_v=zeros(Nom,1);

    ii=0;
    for om = om_v
        ii = ii+1;
   

        % -------------------------------------------- 
        % blocked force from reciprocal approach DFR
           %   (3a.1)  get wavenumbers from dispersion relation 
                    [k0,k] = cal_disproots(para.d,om);
        
           %   (3a.2)  energy to modal density ratio 
                    E2n = pi*para.rho_f*para.g^2/(2*om*k0)*(Sxx(ii)*dom)*(tanh(k0*para.d)+k0*para.d*sech(k0*para.d)^2);%*dom
                    E2n_v(ii) = E2n;

           %   (3a.3) potential damping         
                    Cp = [B11n(ii) B15n(ii);B15n(ii) B55n(ii)];  % 2 by 2 matrix       
                    Cp = Tc.'*Cp*Tc;
           %  (3a.4) block force cross spectrum (DFR)
                    Sfbfb = 4*E2n/pi/om*(om*Cp);
                
                    Sff_om(ii,1)=Sfbfb(1,1); % put the matrix into vectors 
                    Sff_om(ii,2)=Sfbfb(1,2);
                    Sff_om(ii,3)=Sfbfb(2,2);
        
        % -------------------------------------------- 
        % cross spectrum response from reciprocal approach (DFR)
         
            % (3b.1) total damping (add viscous damping)
            Ctot = Cp;
            Ctot(1,1) = Ctot(1,1) + Cv(1,ii);
            Ctot(1,2) = Ctot(1,2) + Cv(2,ii);
            Ctot(2,1) = Ctot(1,2);
            Ctot(2,2) = Ctot(2,2) + Cv(3,ii);
        
            % (3b.2) Added mass         
            Ma = [A11n(ii) A15n(ii);A15n(ii) A55n(ii)];  % 2 by 2 matrix       
            Ma = Tc.'*Ma*Tc;
        
            % (3b.3) Transfer function 
            Hinv = -om^2*(Ms+Ma)+1i*om*Ctot+K;
            H = inv(Hinv);
         
            % (3b.4) structural response cross spectrum (DFR)     
            Szz = H*Sfbfb*H'; % diffused response
                
            Szz_om(ii,1) = abs(Szz(1,1));
            Szz_om(ii,2) = abs(Szz(1,2));
            Szz_om(ii,3) = abs(Szz(2,2));


        % -------------------------------------------- 
        % (3c.1) structural response cross spectrum (Morison's)       
        SyyQ = Hq(:,ii)*(Sxx(ii)*dom)*Hq(:,ii)';
        
        SyyQ_om(ii,1) = abs(SyyQ(1,1));
        SyyQ_om(ii,2) = abs(SyyQ(1,2));
        SyyQ_om(ii,3) = abs(SyyQ(2,2));
        
        % (3c.2) wave forces (Morison's)    
        SffQ = 1/2*Fv(:,ii)*(Sxx(ii)*dom)*Fv(:,ii)'; % ensemble averaged, so divide by 2
        
        SffQ_om(ii,1) = abs(SffQ(1,1));
        SffQ_om(ii,2) = abs(SffQ(1,2));
        SffQ_om(ii,3) = abs(SffQ(2,2));
        
        % (3c.3) separate wave forces into inertia (I) and drag (D) parts (Morison's) 
        SffQ_I = 1/2*FI_v(:,ii)*(Sxx(ii)*dom)*FI_v(:,ii)';
        
        SffQ_I_om(ii,1) = abs(SffQ_I(1,1));
        SffQ_I_om(ii,2) = abs(SffQ_I(1,2));
        SffQ_I_om(ii,3) = abs(SffQ_I(2,2));
        
        SffQ_D = 1/2*FD_v(:,ii)*(Sxx(ii)*dom)*FD_v(:,ii)';
        
        SffQ_D_om(ii,1) = abs(SffQ_D(1,1));
        SffQ_D_om(ii,2) = abs(SffQ_D(1,2));
        SffQ_D_om(ii,3) = abs(SffQ_D(2,2));                    
    end    



%%
% results display

exportfig = 1;

% comparison of forces 
    fig1 = figure;
    
    subplot(121)
    plot(om_v,Sff_om(:,1),'ko-',...
        om_v,SffQ_om(:,1),'r',...
        om_v,SffQ_I_om(:,1),'b.',...
        om_v,SffQ_D_om(:,1),'b--')

    ylim([-1e9 0.8e10])
    xlabel('Frequency [$$rad/s$$]','Interpreter','latex')
    ylabel('Blocked Force [$$N^2$$s/rad]','Interpreter','latex')

    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex','FontSize',16)

    
    subplot(122)
    plot(om_v,Sff_om(:,3),'ko-',...
         om_v,SffQ_om(:,3),'r',...
        om_v,SffQ_I_om(:,3),'b.',...
        om_v,SffQ_D_om(:,3),'b--');
    ylim([-1e13 6e13])

    legend('DFR','Morison''s (Total)','Morison''s (I)','Morison''s (D)',...
        'location', 'NorthEast','Interpreter','latex')
    xlabel('Frequency [$$rad/s$$]','Interpreter','latex')
    ylabel('Blocked Moment [$$N^2m^2s/rad$$]','Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    
    figuresize(30, 12, 'centimeters');
    movegui(fig1, [50 40]);
    set(gcf, 'Color', 'w');


    if exportfig==1
        export_fig FroceComp.png;
    end
    %%    
    % comparison of responese 

    factor_spread=(2*pi*2/pi); % normalising factor to take account of 
    %peak amplitdue from cosine squared spreading sea

    fig1 = figure;
    subplot(121)
    plot(om_v,SyyQ_om(:,1),'r--',...
        om_v,Szz_om(:,1)*factor_spread,'ko-')

    ylim([0 1.5e-4])
    
    xlabel('Frequency [$$rad/s$$]','Interpreter','latex')
    ylabel('Displacement [$$m^2s/rad$$]','Interpreter','latex')
    
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    
    subplot(122)
    plot(om_v,SyyQ_om(:,3),'r--',...
        om_v,Szz_om(:,3)*factor_spread,'ko-')
    ylim([0 5e-7])
    
    legend('Uni-D (Morison''s )','Spreading UBound (DFR)',...
        'location', 'NorthEast','Interpreter','latex')
    xlabel('Frequency [$$rad/s$$]','Interpreter','latex')
    ylabel('Rotation [$$rad^2s/rad$$]','Interpreter','latex')
    
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    
    
    
    
    figuresize(30, 12, 'centimeters');
     movegui(fig1, [50 40])
    set(gcf, 'Color', 'w');

    if exportfig==1
        export_fig ResComp.png;
    end