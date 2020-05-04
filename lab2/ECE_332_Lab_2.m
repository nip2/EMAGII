%% ECE 332 Lab 2
%-------------------------------------------------------------------------%
%--------------------------Lab 2: Waveguides------------------------------%
%-------------------------ver. 1, 03/28/2020------------------------------%
%-------ECE332 Lab, Spring 2020: Special remote education edition---------%

% Lab student co-authors: Nicholas Porter, Regan Garner, Grace Semerjian
% Date: 5/3/2020

%% Clear out junk
clear, clc, close all

%% Part a) Use MATLAB to make the following plots of |~Ez|

% These plots are all of a slice across the waveguide at a fixed z, the 
% Front view below. You can pick any a and Eo, but let b = a/2.  We’re not 
% worried about the numbers, just a clear representation of the field patterns. 

%(i) a 1-d plot of |Ez| vs. x at y = b/2

syms x;
epsilon = 8.854*(10^-12); % electric permittivity / dielectric constant 
mu = 4*pi*(10^7);         % magnetic permeability/ magnetic constant 

a = 1;                    % arbitrary width
b = a/2;                  % length
m = 1;                    % m and n are integers that indicate the mode of 
n = 1;                    % propagation of the TM wave.
y = b/2;                  % fixed value for part (i), looking at 1/2 height of waveguide
Eo = 1;                   % arbitrary strength of the electric field 

Ez_mag(x) = Eo*sin((m*pi*x)/a)*sin((n*pi*y)/b); 
figure(1)
fplot(Ez_mag(x),[0 a])
xlabel('x [m]')
ylabel('|~Ez|')
title('Front View of |Ez| @ y = b/2')
grid on

% (ii) a 2-d contour plot of |Ez| in the xy plane

x = linspace(0,a,100);
y = linspace(0,b,100);
[X,Y] = meshgrid(x,y);
Ez_magxy= Eo*sin((m*pi.*X)/a).*sin((n*pi.*Y)/b); 
figure(2)
contour(X,Y,Ez_magxy)
xlabel('x (0 to a) [m]')
ylabel('y (0 to b) [m]')
title('Front View of |~Ez|')
grid on

% (iii) a 3-d surface plot of |Ez| in the xy plane

x = linspace(-a,a,100);
y = linspace(0,b,100);
[X,Y] = meshgrid(x,y);
Ez_magxy= Eo*sin((m*pi.*X)/a).*sin((n*pi.*Y)/b); 
figure(3)
surface(X,Y,Ez_magxy)
xlabel('x (-a to a) [m]')
ylabel('y (0 to b) [m]')
title('Front View of |~Ez|')
grid on
colorbar
view(3)

%Comments on the 3 figures
%The 1-D, 2-D, and 3-D figures all seem to agree in the plot for |~Ez|. We
%can judge this by the outlines, values, and colors. If we compare the 1-D
%plot versus the 3-D plot we see that the cross section checks out with the
%surface plot with the x values vs the |~Ez| values @ y = b/2. Note that
%for the 3-D plot x was from [-a,a] but our cross-section is from [0,a]. If
%we then compare the 2-D contour plot with the 3-D surface plot we see that
%again, the placements line up as do their colors that denote their values.
%Again, note that for the 3-D plot x was from [-a,a] but our contour plot
%is from [0,a]. Lastly, as a brief sanity check, we should be expecting a 
%sinusoidal function as our plot for |~Ez|, and that is what we see here.


%% Part b) Create three additional plots. 
% One should be a different mode, one a different a and b ratio, and one 
% should be |~Hz|.  Exactly what you do here – which mode, what ratio, what 
% type of plot – is up to you. Experiment a bit, and turn in whatever you 
% think is interesting and informative.

% Plot of a different mode. Let m = 8 & n = 24

syms x;
epsilon = 8.854*(10^-12); % electric permittivity / dielectric constant 
mu = 4*pi*(10^7);         % magnetic permeability/ magnetic constant 

a = 1;                    % arbitrary width
b = a/2;                  % length
m = 8;                    % m and n are integers that indicate the mode of 
n = 24;                   % propagation of the TM wave.
Eo = 1;                   % arbitrary strength of the electric field 

x = linspace(-a,a,100);
y = linspace(0,b,100);
[X,Y] = meshgrid(x,y);
Ez_magxy= Eo*sin((m*pi.*X)/a).*sin((n*pi.*Y)/b); 
figure(4)
surface(X,Y,Ez_magxy)    % What abomination should come out?
xlabel('x (-a to a) [m]')
ylabel('y (0 to b) [m]')
title('Front View of |~Ez|')
grid on
colorbar
view(3)

% Plot of a different a & b ratio 

a = 50;                   % arbitrary width
b = 1;                    % arbitrary length
                    
m = 1;                    % m and n are integers that indicate the mode of 
n = 1;                    % propagation of the TM wave.
Eo = 1;                   % arbitrary strength of the electric field 

x = linspace(-a,a,100);
y = linspace(-b,b,100);
[X,Y] = meshgrid(x,y);
Ez_magxy= Eo*sin((m*pi.*X)/a).*sin((n*pi.*Y)/b); 
figure(5)
surface(X,Y,Ez_magxy)    % What abomination should come out?
xlabel('x (-a to a) [m]')
ylabel('y (-b to b) [m]')
title('Front View of |~Ez|')
grid on
colorbar
view(3)

% Plot of |~Hz| for a TE wave 

a = 1;                    % arbitrary width
b = a/2;                  % length
m = 1;                    % m and n are integers that indicate the mode of 
n = 1;                    % propagation of the TM wave.
Ho = 1;                   % arbitrary strength of the magnetic field 

x = linspace(-a,a,100);
y = linspace(-b,b,100);
[X,Y] = meshgrid(x,y);
Hz_magxy = Ho*cos((m*pi.*X)/a).*cos((n*pi.*Y)/b); 
figure(6)
surface(X,Y,Hz_magxy)    
xlabel('x (-a to a) [m]')
ylabel('y (-b to b) [m]')
title('Front View of |~Hz|')
grid on
colorbar
view(3)

%Comments on the 3 figures
%For the first figure in Part b) we're playing around with changing the
%mode (we let m = 8 & n = 24) while leaving all of the other paramters from
%part a the same for the |~Ez| phasor component. We went with a surface
%plot to gain a more bird's eye view of what's going on with the TM wave.
%As we can see, it gets a whole lot more cluttered with higher frequency
%sinusoidal waves, but remain within the same amplitude as we saw in Pt a).
%This makes sense as the arguments for |~Ez| sin()'s are much larger than 
%before. For the second plot we messed around with changing the a and b 
%ratio while keeping everything from Part a) the same. Judging from the 
%provided surface plot (and the other plots we had when playing around with 
%other values), it appears that even when a and b vary, the overall shape 
%and amplitude that we saw in Pt. a) does not change. Note that for this
%plot x was from [-a,a] and y is from [-b,b]. It's merely scaled to
%the larger (or smaller a and b) values over space. Again, this is expected
%if we look in the arguments for the sin()'s in |~Ez|. For the third figure
%we plotted |~Hz| for a TE wave (as |~Hz| = 0 for a TM wave) with the same
%parameters as we had in Pt. a). We should once more be expecting a
%sinusoidal function with a max amplitude of 1 (cnosen Ho) for the output, 
%which is indeed what we see on our surface plot.


%% Part c) Write a code for one of the plots where the user has to choose some options.  
% This involves a prompt to the screen where the user would enter some 
% information, e.g., which mode, or a value for a, or which type of plot.
% Show a couple of example outputs.

%Comments for Part c)
%We decided to tie all of what we learned about waveguides together in one
%go by asking the user for which wave type (TE or TM),a,b,mode,carrier
%frequency, field strength, phasor vector magnitude component, & what type 
%of plot (1-D cross-sectional (spits out both side and front views @  
%y = b/2 and x = a/2), 2-D contour, or 3-D surface). When going through the
%questions, the user is prompted to put in feasible values such as a
%carrier frequency that is => the cutoff frequency fmn, or the modes that
%work given the selected a and b values.

%All of what's set up below (parameters, directions, other conventions, etc.)
%is wrt the diagram that was provided for this lab.

% TM, or TE?
wave_type = input('Which wave type, TM or TE?: \n', 's');

while ~(strcmp(wave_type, 'TM') | strcmp(wave_type, 'TE')) 
    disp('The question was not answered properly')
    wave_type = input('Please enter a wavetype, TM or TE:\n ', 's');
end

% What a & b?
a = input('What is the preferred width (x direction) of the waveguide?: \n');
b = input('What is the preferred height (y direction) of the waveguide?: \n');

while a == 0 | b == 0
        disp('2 dimensions are needed')
        if a == 0
        a = input('Please enter a value that is not equal to 0: \n');
        end
        if b == 0
            b = input('Please enter a value that is not equal to 0: \n');
        end
end

abigb = a > b;

% What mode?
disp('What mode is preferred?')
mode_m = input('Please enter a positive integer value for m that is greater than or equal to 1: \n');
    
while ~(floor(mode_m)==ceil(mode_m)) | mode_m < 1 & abigb ~= 1 
        disp('This is not a valid value for m')
        mode_m = input('Please enter non-negative integer value for m (if a > b, m must be greater than or equal to 1): \n');
end

mode_n = input('Please enter a positive integer value for n (must be => 1 if TM): \n');
    
while ~(floor(mode_n)==ceil(mode_n)) | mode_n < 0 | (mode_n < 1 & wave_type == 'TM')
        disp('This is not a valid value for n')
        mode_n = input('Please enter a positive integer value for n (must be => 1 if TM): \n');
end

%Trying to figure out the cutoff frequency fmn 
epsilon = 8.854*(10^-12); % electric permittivity / dielectric constant 
mu = 4*pi*(10^7);         % magnetic permeability/ magnetic constant
m = mode_m;
n = mode_n;
c = 3e8;  % We're assuming that the waveguide is full of vacuum, so upo = c
fmn = (c/2)*sqrt((m/a)^2 + (n/b)^2);

%What frequency?
disp('What frequency should the wave be working at?')
f = input(['Please enter a number in Hz that is greater than or equal to the cutoff frequency ' num2str(fmn) ': \n']);

while ~(f >= fmn)
        disp('This is not a valid value for the frequency')
        f = input(['Please enter a value that is greater than or equal to the cutoff frequency ' num2str(fmn) ': \n']);
end

%What field strength?
strength = input('What field strength should the wave be able to max out at? Please enter a positive number: \n');

while ~(strength >= 0)
        disp('This is not a valid value for the field strength')
        strength = input('Please enter a positive value for the field strength: \n');
end

% What phasor vector magnitude component?
disp('Which phasor vector magnitude component would you like to plot?')
component = input('Enter 1 for |~Ex|, 2 for |~Ey|, 3 for |~Ez|, 4 for |~Hx|, 5 for |~Hy|, 6 for |~Hz|: \n');

while ~(component == 1 | component == 2 | component == 3 | component == 4 | component == 5 | component == 6)
        disp('Not a valid choice')
        component = input('Please enter a 1 for |~Ex|, 2 for |~Ey|, 3 for |~Ez|, 4 for |~Hx|, 5 for |~Hy|, 6 for |~Hz|: \n');
end

% What plot?
disp('What kind of plot is preferred?')
plot_type = input('Enter 1 for a 1-D cross-section plot, 2 for a 2-D contour plot, 3 for a 3-D surface plot: \n');

while ~(plot_type == 1 | plot_type == 2 | plot_type == 3)
        disp('Not a valid choice')
        plot_type = input('Please enter a 1 for 1-D cross-section plot, 2 for a 2-D contour plot, 3 for a 3-D surface plot: \n');
end

%Setting up some terms for calculating the desired phasor magnitude comp.'s
w = 2*pi*f;
mpi = m*pi/a;
npi = n*pi/b;
beta = sqrt(w^2*mu*epsilon - mpi^2 - npi^2);
k = w*sqrt(mu*epsilon); %Note: all of this is assumed to be in vacuum
kc_sq = k^2 - beta^2;
Eo = strength;
Ho = strength;

bkc = beta/kc_sq;  %got rid of i's
wekc = w*epsilon/kc_sq;
wmukc = w*mu/kc_sq;

if plot_type == 1 %Plotting 1-D cross-sections
    
    if wave_type == 'TM'
        
        if component == 1 %If |~Ex| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Ex_mag(x) = -bkc*mpi*Eo*cos(mpi*x)*sin(npi*y);
        figure(7)
        fplot(Ex_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Ex|')
        title('Front View of |~Ex| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Ex_mag(y) = -bkc*mpi*Eo*cos(mpi*x)*sin(npi*y);
        figure(8)
        fplot(Ex_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Ex|')
        title('Front View of |~Ex| @ x = a/2')
        grid on
        end
        
        if component == 2 %If |~Ey| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Ey_mag(x) = -bkc*npi*Eo*sin(mpi*x)*cos(npi*y);
        figure(7)
        fplot(Ey_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Ey|')
        title('Front View of |~Ey| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Ey_mag(y) = -bkc*npi*Eo*sin(mpi*x)*cos(npi*y);
        figure(8)
        fplot(Ey_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Ey|')
        title('Front View of |~Ey| @ x = a/2')
        grid on
        end
        
        if component == 3 %If |~Ez| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Ez_mag(x) = Eo*sin(mpi*x)*sin(npi*y); 
        figure(7)
        fplot(Ez_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Ez|')
        title('Front View of |Ez| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Ez_mag(y) = Eo*sin(mpi*x)*sin(npi*y);
        figure(8)
        fplot(Ez_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Ez|')
        title('Front View of |~Ez| @ x = a/2')
        grid on
        end
        
        if component == 4 %If |~Hx| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Hx_mag(x) = wekc*npi*Eo*sin(mpi*x)*cos(npi*y);
        figure(7)
        fplot(Hx_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Hx|')
        title('Front View of |~Hx| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Hx_mag(y) = wekc*npi*Eo*sin(mpi*x)*cos(npi*y);
        figure(8)
        fplot(Hx_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Hx|')
        title('Front View of |~Hx| @ x = a/2')
        grid on
        end
        
        if component == 5 %If |~Hy| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Hy_mag(x) = -wekc*mpi*Eo*cos(mpi*x)*sin(npi*y);
        figure(7)
        fplot(Hy_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Hy|')
        title('Front View of |~Hy| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Hy_mag(y) = wekc*npi*Eo*sin(mpi*x)*cos(npi*y);
        figure(8)
        fplot(Hy_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Hy|')
        title('Front View of |~Hy| @ x = a/2')
        grid on
        end
        
        if component == 6 %If |~Hz| is chosen
        % Hz_magxy = 0;
        disp('This is a TM wave, so there is no |~Hz| component to plot')
        end           
    end
      
    if wave_type == 'TE' % For a TE wave      
        
        if component == 1 %If |~Ex| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Ex_mag(x) = wmukc*npi*Ho*cos(mpi*x)*sin(npi*y);
        figure(7)
        fplot(Ex_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Ex|')
        title('Front View of |~Ex| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Ex_mag(y) = wmukc*npi*Ho*cos(mpi*x)*sin(npi*y);
        figure(8)
        fplot(Ex_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Ex|')
        title('Front View of |~Ex| @ x = a/2')
        grid on
        end
        
        if component == 2 %If |~Ey| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Ey_mag(x) = -wmukc*mpi*Ho*sin(mpi*x)*cos(npi*y);
        figure(7)
        fplot(Ey_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Ey|')
        title('Front View of |~Ey| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Ey_mag(y) = -wmukc*mpi*Ho*sin(mpi*x)*cos(npi*y);
        figure(8)
        fplot(Ey_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Ey|')
        title('Front View of |~Ey| @ x = a/2')
        grid on
        end
        
        if component == 3 %If |~Ez| is chosen  
        disp('This is a TE wave, so there is no |~Ez| component to plot')
        end
        
        if component == 4 %If |~Hx| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Hx_mag(x) = bkc*mpi*Ho*sin(mpi*x)*cos(npi*y);
        figure(7)
        fplot(Hx_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Hx|')
        title('Front View of |~Hx| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Hx_mag(y) = bkc*mpi*Ho*sin(mpi*x)*cos(npi*y);
        figure(8)
        fplot(Hx_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Hx|')
        title('Front View of |~Hx| @ x = a/2')
        grid on
        end
        
        if component == 5 %If |~Hy| is chosen  
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Hy_mag(x) = bkc*npi*Ho*cos(mpi*x)*sin(npi*y);
        figure(7)
        fplot(Hy_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Hy|')
        title('Front View of |~Hy| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Hy_mag(y) = bkc*npi*Ho*cos(mpi*x)*sin(npi*y);
        figure(8)
        fplot(Hy_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Hy|')
        title('Front View of |~Hy| @ x = a/2')
        grid on
        end
        
        if component == 6 %If |~Hz| is chosen
        %Front view cross-section @ y = b/2
        syms x;
        y = b/2;
        Hz_mag(x) = Ho*cos(mpi*x)*cos(npi*y); 
        figure(7)
        fplot(Hz_mag(x),[0 a])
        xlabel('x [m]')
        ylabel('|~Hz|')
        title('Front View of |~Hz| @ y = b/2')
        grid on
        
        %Side view cross-section @ x = a/2
        syms y;
        x = a/2;
        Hz_mag(y) = Ho*cos(mpi*x)*cos(npi*y); 
        figure(8)
        fplot(Hz_mag(y),[0 b])
        xlabel('y [m]')
        ylabel('|~Hz|')
        title('Front View of |~Hz| @ x = a/2')
        grid on
        end
    end 
end

if plot_type == 2 %Plotting 2-D contour plots
    
    if wave_type == 'TM' % For a TM wave
        
        if component == 1 %If |~Ex| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Ex_magxy = -bkc*mpi*Eo*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        contour(X,Y,Ex_magxy)    
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Ex|')
        grid on
        end
        
        if component == 2  %If |~Ey| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Ey_magxy = -bkc*npi*Eo*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        contour(X,Y,Ey_magxy)    
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Ey|')
        grid on
        end
        
        if component == 3 %If |~Ez| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Ez_magxy = Eo*sin(mpi.*X).*sin(npi.*Y);
        figure(7)
        contour(X,Y,Ez_magxy)    
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Ez|')
        grid on
        end

        if component == 4 %If |~Hx| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Hx_magxy = wekc*npi*Eo*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        contour(X,Y,Hx_magxy)   
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Hx|')
        grid on
        end

        if component == 5 %If |~Hy| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Hy_magxy = -wekc*mpi*Eo*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        contour(X,Y,Hy_magxy)   
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Hy|')
        grid on
        end

        if component == 6 %If |~Hz| is chosen
        % Hz_magxy = 0;
        disp('This is a TM wave, so there is no |~Hz| component to plot')
        end
    end
    
    if wave_type == 'TE' % For a TE wave
      
        if component == 1 %If |~Ex| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Ex_magxy = wmukc*npi*Ho*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        contour(X,Y,Ex_magxy)    
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Ex|')
        grid on
        end

        if component == 2  %If |~Ey| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Ey_magxy = -wmukc*mpi*Ho*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        contour(X,Y,Ey_magxy)    
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Ey|')
        grid on
        end
        
        if component == 3 %If |~Ez| is chosen
        % Ez_magxy = 0;    
        disp('This is a TE wave, so there is no |~Ez| component to plot')
        end

        if component == 4 %If |~Hx| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Hx_magxy = bkc*mpi*Ho*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        contour(X,Y,Hx_magxy)   
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Hx|')
        grid on
        end

        if component == 5 %If |~Hy| is chosen
        x = linspace(0,a,100);
        y = linspace(0,b,100);
        [X,Y] = meshgrid(x,y);
        Hy_magxy = bkc*npi*Ho*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        contour(X,Y,Hy_magxy)   
        xlabel('x (0 to a) [m]')
        ylabel('y (0 to b) [m]')
        title('Front View of |~Hy|')
        grid on
        end

        if component == 6 %If |~Hz| is chosen
        Hz_magxy= Ho*cos((m*pi.*X)/a).*cos((n*pi.*Y)/b);
        end        
    end   
end

if plot_type == 3 %Plotting 3-D surface plots
    
    if wave_type == 'TM'
        
        if component == 1
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Ex_magxy = -bkc*mpi*Eo*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        surface(X,Y,Ex_magxy)    % What abomination should come out?
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Ex|')
        grid on
        colorbar
        view(3)
        end

        if component == 2  %If |~Ey| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Ey_magxy = -bkc*npi*Eo*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        surface(X,Y,Ey_magxy)    
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Ey|')
        grid on
        colorbar
        view(3)
        end
        
        if component == 3 %If |~Ez| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Ez_magxy = Eo*sin(m*pi.*X).*sin(n*pi.*Y);
        figure(7)
        surface(X,Y,Ez_magxy)    
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Ez|')
        grid on
        colorbar
        view(3)
        end

        if component == 4 %If |~Hx| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Hx_magxy = wekc*npi*Eo*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        surface(X,Y,Hx_magxy)   
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Hx|')
        grid on
        colorbar
        view(3)
        end

        if component == 5 %If |~Hy| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Hy_magxy = -wekc*mpi*Eo*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        surface(X,Y,Hy_magxy)   
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Hy|')
        grid on
        colorbar
        view(3)
        end

        if component == 6 %If |~Hz| is chosen
        % Hz_magxy = 0;
        disp('This is a TM wave, so there is no |~Hz| component to plot')
        end
    end
    
    if wave_type == 'TE' % For a TE wave

        if component == 1
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Ex_magxy = wmukc*npi*Ho*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        surface(X,Y,Ex_magxy)    % What abomination should come out?
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Ex|')
        grid on
        colorbar
        view(3)
        end

        if component == 2  %If |~Ey| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Ey_magxy = -wmukc*mpi*Ho*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        surface(X,Y,Ey_magxy)    
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Ey|')
        grid on
        colorbar
        view(3)
        end
        
        if component == 3 %If |~Ez| is chosen
        % Ez_magxy = 0;    
        disp('This is a TE wave, so there is no |~Ez| component to plot')
        end

        if component == 4 %If |~Hx| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Hx_magxy = bkc*mpi*Ho*sin(mpi.*X).*cos(npi.*Y);
        figure(7)
        surface(X,Y,Hx_magxy)   
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Hx|')
        grid on
        colorbar
        view(3)
        end

        if component == 5 %If |~Hy| is chosen
        x = linspace(-a,a,100);
        y = linspace(-b,b,100);
        [X,Y] = meshgrid(x,y);
        Hy_magxy = bkc*npi*Ho*cos(mpi.*X).*sin(npi.*Y);
        figure(7)
        surface(X,Y,Hy_magxy)   
        xlabel('x (-a to a) [m]')
        ylabel('y (-b to b) [m]')
        title('Front View of |~Hy|')
        grid on
        colorbar
        view(3)
        end

        if component == 6 %If |~Hz| is chosen
        Hz_magxy= Ho*cos((m*pi.*X)/a).*cos((n*pi.*Y)/b);
      
        end

    end  
end

    
% %TM wave's magnitude phasor vector components
% Ex_magxy = -bkc*mpi*Eo*cos(mpi.*x).*sin(npi.*y);
% Ey_magxy = -bkc*npi*Eo*sin(mpi.*x).*cos(npi.*y);
% Hx_magxy = wekc*npi*Eo*sin(mpi.*x).*cos(npi.*y);
% Hy_magxy = -wekc*mpi*Eo*cos(mpi.*x).*sin(npi.*y);
% Ez_magxy = Eo*sin(mpi.*X).*sin(npi.*Y); 
% % Hz_magxy = 0;
% 
% %TE wave's magnitude phasor vector components
% Ex_magxy = wmukc*npi*Ho*cos(mpi.*x).*sin(npi.*y);
% Ey_magxy = -wmukc*mpi*Ho*sin(mpi.*x).*cos(npi.*y);
% Hx_magxy = bkc*mpi*Ho*sin(mpi.*x).*cos(n*pi.*y);
% Hy_magxy = bkc*npi*Ho*cos(mpi.*x).*sin(n*pi.*y);
% Hz_magxy= Ho*cos(mpi.*X).*cos(n*pi.*Y); 
% % Ez_magxy = 0;
