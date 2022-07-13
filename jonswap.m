function [ S, Amp, Phase ] = jonswap( Ohm, Hs, Tp)
%JONSWAP - Calculates the wave spectrum values for a JONSWAP spectrum

% Ohm - wave frequency [rad/s]
% Hs - significant wave height
% Tp - wave period associated with the peak of the spectrum 

% 21/03/2019 @ JDG-01, Cambridge  [J Yang] 

wp = 2*pi/Tp; % peak frequency 
Gamma = 1;  % peakedness factor 


N=length(Ohm);
S=zeros(N,1);
 for ii = 1:N
     if Ohm(ii)<wp
         Sigma = 0.07;
     else
         Sigma = 0.09;
     end
     A = exp(-((Ohm(ii)/wp-1)/(Sigma*sqrt(2)))^2);
     
     if Ohm(ii)==0
         S(ii)=0;
     else
         S(ii) = 320*Hs^2*Ohm(ii)^-5/Tp^4*exp(-1950*Ohm(ii)^-4/Tp^4)*Gamma^A;
%          S(ii) = 0.0081*9.8^2*Ohm(ii)^-5*exp(-1.25*Ohm(ii)^-4*(2*pi/Tp)^4)*Gamma^A;
     end
 end

 % Determine the frequency step from the frequency vector. Note that the
 % highest frequency step is extrapolated.
 domg = zeros(N,1);
 domg(1:end-1) = diff( Ohm );
 domg(end) = domg(end-1);

 % Determine the amplitudes from the spectral values
 Amp = sqrt( 2 * S .* domg );

 % Random phases
 Phase = rand(length(Ohm),1)*2*pi;

end