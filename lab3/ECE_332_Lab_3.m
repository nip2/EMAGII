%% ECE 332 Lab 3
%-------------------------------------------------------------------------%
%--------------------------Lab 3: Antennas--------------------------------%
%-------ECE332 Lab, Spring 2020: Special remote education edition---------%
%-------------------------------------------------------------------------%

% Lab student co-authors: Regan Garner, Grace Semerjian
% Date: 5/17/2020

%% Clear out junk
clear, clc, close all

%% Part 1) Plot ~I(z) vs z for a dipole antenna of length L = 1.25*lambda 
%oriented in the z-direction. Assume z = 0 is the middle of the antenna, so
%your plot goes from z = -L/2 to L/2.  Let Io = 1.

% Setting up I(z)

syms z

lambda = 1;                             %Let lambda = 1 here
l = 1.25*lambda;                        %L = 1.25*lambda is given
Io = 1;                                 %Io = 1 is given
k = (2*pi)/lambda;

Itophalf = Io*sin(k*(l/2 - z));         %~I(z) for 0 <= z <= l/2
Ibothalf = Io*sin(k*(l/2 + z));         %~I(z) for -l/2 <= z < 0

figure(1)                               %Plotting ~I(z)'s piece-wise func.
fplot(Itophalf,z,[0 l/2])              
hold on
fplot(Ibothalf,z,[-l/2 0])
grid on
title('~I(z) for L = 1.25*lambda Dipole Antenna')
xlabel('~I(z) [A]')
ylabel('z [m]')

%Comments for Part 1)
%As we can see with the figure created, ~I(z) resembles the intermediate
%plot of L = 1*lambda and L = 1.5*lambda. This is expected since ours is L
%= 1.25* lambda. As another check we see that ~I(z) does indeed go to zero
%at z = -/+ L/2


%% Part 2) Generate a polar plot of the normalized radiation intensity 
%F(theta) for the antenna in step 1. You will need the avg. power density
%S(theta). This is given in equation 9.56 in the 7th edition of Ulaby. 
%Since we are looking at an arbitrary fixed distance R, we can simply let 
%the constant So = 15Io^2/piR^2 = 1 (it cancels out in the next step anyway)
%Use MATLAB to find the maximum of S(theta), and calculate F(theta) = 
%S(theta)/Smax. Note that Smax ~= So for any L. Also be careful that you are 
%using theta in radians. As a check, your plot should look something like 
%Fig. 9-17 (b), but with a few more small lobes.  If you want, try running
%your program with L = 0.5*lambda, lambda, & 1.5*lamda. You should generate 
%Fig. 9-17 exactly.

lambda = 1;                              %Arbitrarily let lambda = 1
theta1 = linspace(0,pi,100);             %Would like to plot against both
theta2 = linspace(0,-pi,100);            %sides of the z axis
Io = 1;                                  %Io = 1 is given
R = 1;                                   %Arbitrarily let R = 1
So = (15*Io^2)/(pi*R^2);                 %So = 15Io^2/piR^2 = 1 is given

%Just trying to match Figure 9-17's plots
%l = lambda/2;
%pit = (pi*l)/lambda;
% Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
% Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
% figure(2)
% polar(theta1,Savg1);
% hold on
% polar(theta2,Savg2);
% title('(Not normalized) Radiation pattern for L = lambda/2')
% 
% l = lambda;
% pit = (pi*l)/lambda;
% Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
% Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
% figure(3)
% polar(theta1,Savg1);
% hold on
% polar(theta2,Savg2);
% title('(Not normalized) Radiation pattern for L = lambda')

%This is the plot that is of interest (from Part 1)
l = 1.25*lambda;
pit = (pi*l)/lambda; %Just a constant that I didn't want to rewrite ea time
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
%The 'pi/2 - theta#' terms seen in the cos() and sin() arguments are merely
%to match the orientation as seen in the textbook
figure(4)
polar(theta1,Savg1);
hold on
polar(theta2,Savg2);
title('(Not normalized) Radiation pattern for L = 1.25*lambda')

% l = 1.5*lambda;
% pit = (pi*l)/lambda;
% Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
% Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
% figure(5)
% polar(theta1,Savg1);
% hold on
% polar(theta2,Savg2);
% title('(Not normalized) Radiation pattern for L = 1.5*lambda')

%Plots check out, carrying on with finding F(theta)
Savg1_max = max(Savg1);   %Finding Smax for both sides of the axis, should 
Savg2_max = max(Savg2);   %be the same for both sides
Fone = Savg1/Savg1_max;   %Using Smax to find F(theta)
Ftwo = Savg2/Savg2_max;
figure(6)        
polar(theta1,Fone);       %Plotting the normalized avg. power density for
hold on                   %L = 1.25*lambda
polar(theta2,Ftwo);
title('(Normalized) Radiation pattern for L = 1.25*lambda')

%To be used for Part 3 calculation
Savg1 = So*((cos(pit.*cos(theta1)) - cos(pit))./sin(theta1)).^2;
Savg1_max = max(Savg1);
%This calculation was rerun again because of the 'pi/2 - theta#' terms from
%before, they were not computationally correct, only aesthecially correct

%Comments for Part 2)
%Again, our F(theta) plot for the L = 1.25*lambda plot should look like the
%intermediate plot of the 1*lambda and 1.5*lambda plots from Fig. 9-17.
%Ours does indeed look like it's supposed to with the small 4 leaf clover 
%lobes centered around the origin and the main lobes (peak at 1) centered
%around theta = 90 degrees. Note that the commented out plots were the L =
%0.5*lambda, 1*lambda, and 1.5*lambda plots as something to check against.
%We also didn't realize that for Part 4) we'd be calculating S(theta) and
%just did it prematurely here. 


%% Part 3) The pattern solid angle OmegaP and directivity D are measures of 
%how narrow an antenna’s radiation pattern is. They are related by D = 
%4pi/OmegaP. An isotropic antenna has OmegaP = 4pi and D = 1. The narrower 
%the radiation pattern, the smaller OmegaP the larger D. Find OmegaP and D 
%by numerically integrating F(theta) for the 1.25*lambda antenna in the 
%steps above. The integration can be approximated by OmegaP approximately = 
%2pi*sum(Fi(theta))*delta(theta). You’ll need to define an increment of F 
%for each increment of theta. The smaller the increment, the more accurate 
%the approximation. Answers will vary a bit, but you should get OmegaP = 
%3.8. Alternately, you can use an integration routine in MATLAB if you prefer. 

%Going with integration method
%OmegaP = int(int(F(theta,phi),dOm));

syms theta                   

lambda = 1;          %Arbitrarily let lambda = 1
l = 1.25*lambda;     %Our dipole antenna length of interest
pit = (pi*l)/lambda; %Just a constant that I didn't want to rewrite ea time
Io = 1;              %Io = 1 is given
R = 1;               %Arbitrarily let R = 1
So = (15*Io^2)/(pi*R^2); %So is given 
%Recalculating the various terms for integration
Savg(theta) = So*((cos(pit*cos(theta)) - cos(pit))/sin(theta)).^2;
F(theta) = Savg(theta)/Savg1_max;
OmegaP = vpaintegral(2*pi*(F(theta)*sin(theta)),theta,[0 pi]);
%Utilized vpaintegral() since int() left the integral in it's unsolved form
fprintf('The pattern solid angle OmegaP for the 1.25*lambda dipole antenna is: ')
disp(OmegaP)
D = (4*pi)/OmegaP;
fprintf('The directivity D for the 1.25*lambda dipole antenna is: %4.4f\n',D)

%Comments for Part 3)
%Our L = 1.25*lambda dipole antenna has more directivity than the 0.5*lambda
%antenna as expected (visually it's more elongated at its peak, beta, and
%thus its OmegaP is smaller too)


%% Part 4) Lastly, make a series of polar plots of S(theta) for L from
%0.1*lambda to 2.1*lambda. You can adjust the increments of L as you want, 
%but using increments of 0.1*lambda gives 20 steps and works well. Note we 
%are looking at S, not the normalized F, here so that we can see the 
%different power levels as well as the radiation patterns as the length is 
%varied. The calculation of S(theta) is the same as in step 2 above. At the
%end of the loop for each value of L, insert a pause command. When you run 
%the program, it should display the first plot, then when you hit enter it 
%will show the next one, and so on. You can manually step through each plot
%Put the plots for L = 0.5*lambda, 1*lambda, and 1.5*lambda in your report.

%Same as before
lambda = 1;
theta1 = linspace(0,pi,100);
theta2 = linspace(0,-pi,100);
Io = 1;
R = 1;
So = (15*Io^2)/(pi*R^2);
lens = [0.1:0.1:2.1];       %Each lens array element is to be multiplied
                            %against lambda to make our desired L

for i = 1:length(lens)      %Looking at each specified L polar plot with 
                            %the pause command
inc = lens(i);
l = inc*lambda;
pit = (pi*l)/lambda;
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
figure(7)
polar(theta1,Savg1);
hold on
polar(theta2,Savg2);
title(['(Not normalized) Radiation pattern for L = ' num2str(inc) '*lambda'])
pause
close all

end

%To be included in the report
%L = 0.5*lambda dipole not normalized polar plot
l = lambda/2;
pit = (pi*l)/lambda;
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
figure(8)
polar(theta1,Savg1);
hold on
polar(theta2,Savg2);
title('(Not normalized) Radiation pattern for L = lambda/2')

%L = 1*lambda dipole not normalized polar plot
l = lambda;
pit = (pi*l)/lambda;
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
figure(9)
polar(theta1,Savg1);
hold on
polar(theta2,Savg2);
title('(Not normalized) Radiation pattern for L = lambda')

%L = 1.5*lambda dipole not normalized polar plot
l = 1.5*lambda;
pit = (pi*l)/lambda;
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
figure(10)
polar(theta1,Savg1);
hold on
polar(theta2,Savg2);
title('(Not normalized) Radiation pattern for L = 1.5*lambda')

%Comments for Part 4)
%The not normalized plots should follow the normalized plots as far as
%shape goes, and it appears that our plots follow that. 


%% Extra Credit. Using the command getframe in your program, you should be 
%able to save each plot, and then use the command movie to run them as an 
%animation. I didn’t get this to work with a quick try, but if you explore 
%these commands it should be possible to play back the plots of the 
%radiation pattern as you watch. You will want to use smaller increments of 
%L for a smooth animation. This will be worth up to 10 points extra credit.

%Same as before
lambda = 1;
theta1 = linspace(0,pi,100);
theta2 = linspace(0,-pi,100);
Io = 1;
R = 1;
So = (15*Io^2)/(pi*R^2);
lens = [0.1:0.05:2.1]; %Each lens array element is to be multiplied against 
                       %lambda to make our desired L
F = getframe(gcf);     %Initializing F for whole code to run

for i = 1:length(lens)  %Creating the snapshots for each specified L

inc = lens(i);
l = inc*lambda;
pit = (pi*l)/lambda;
Savg1 = So*((cos(pit.*cos(pi/2 - theta1)) - cos(pit))./sin(pi/2 - theta1)).^2;
Savg2 = So*((cos(pit.*cos(pi/2 - theta2)) - cos(pit))./sin(pi/2 - theta2)).^2;
figure(11)
polarplot(theta1,Savg1);
hold on
polarplot(theta2,Savg2);
hold off                %this little booger line of code made all the diff
rlim([0 30])
title(['(Not normalized) Radiation pattern for L = ' num2str(inc) '*lambda'])
F(i) = getframe(gcf);
 
end

%The animation itself
movie(figure(12),F)

%Comments for Extra Credit
%Man that was cool, still wish I understood the math a bit better to
%correlate to the pictures more easily.

