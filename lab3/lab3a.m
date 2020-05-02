% ECE332: lab 3
clc
close all
clear all
count = 1;

% part 1)
% an antenna of length L, along the z axis with center at the origin
% a current, at microwave range of frequency, through the antenna

c = 2.99792458*10^8;% speed of light
I0 = 1;             % current amplitude
f = 2.4*10^9;       % frequency of both a microwave oven and wifi
lambda = c/f;       % wavelength
k = (2*pi)/lambda;  % wave #
L = 1.25*lambda;    % length of antenna
N = 100;            % # of data points
z = linspace(-L/2,L/2,N);

% piecewise equation for current phasor I~(z)
%{
for i = 1:N
    if z(i) < 0
        I(i) = I0*sin(k*(L/2+z(i)));
    elseif abs(z(i)) > L/2
        disp('error: z out of bounds')
    else
        I(i) = I0*sin(k*(L/2-z(i)));
    end
end
%}
% the above code works, but is bulky. use below instead:
I = I0.*sin(k.*(L/2-abs(z)));

% plot phasor current distribution I(z) vs position on antenna, z
% this plot looks right because the length of the antenna is 
% shorter than 1.5*lambda of the signal frequency, so the 
% full waveform doesn't fit

figure(count)
plot(z,I)
title('$\widetilde{I}$(z) vs z','Interpreter','latex')
xlabel('z, length of antenna (m)','Interpreter','latex')
ylabel('$\widetilde{I}$(z) (A)','Interpreter','latex')
grid minor

count = count + 1;

%%
% part 2)


eta_0 = 120*pi;              % intrinsic impedance
R = 1;                       % distance from the z-axis
theta = linspace(-pi,pi,N);
E_theta = 1i*60*I0*(exp(-1i*k*R)/R).*((cos(k*L/2.*cos(theta))-cos(k*L/2))./sin(theta));
ab_E = (abs(E_theta)).^2;
S = ab_E./(2*eta_0);
S_max = max(S, [], 'all');
F = S./S_max;                % normalized radiation intensity
%alternatively, we should see if this is equivalent 
% F = [(cos((pi*L)\lambda)*cos(theta))-cos((pi*L/lambda)/(sin(theta)))]^2;


figure(count)
plot(theta,F)
title('$F(\theta)(z)$ vs z','Interpreter','latex');
xlabel('$z$, length of antenna $(m)$','Interpreter','latex');
ylabel('$F(\theta) (A)$','Interpreter','latex');
grid minor

count = count + 1;

%%
% part 3)
% The pattern solid angle Omega_p and directivity D are measures of how narrow an antenna’s 
% radiation pattern is.  They are related by D = 4*pi/Omega_p.  An isotropic antenna has 
% Omega_p = 4*pi and D = 1.  The narrower the radiation pattern, the smaller pand the larger D.
% Find pand D by numerically integrating F(theta) for the 1.25*lambda antenna in the steps above.  
% The integration can be approximated by Omega_p≈ 2Fi().  
% You’ll need to define an increment of F for each increment of .  
% The smaller the increment, the more accurate the approximation.  
% Answers will vary a bit, but you should get Omega_p≈ 3.8. 
% Alternately, you can use an integration routine in MATLAB if you prefer.


