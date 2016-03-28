% FLUVIAL

type = 'exner'; % (Multiple options; uses a switch)
    % exner = transport-limited with Exner continuity equation
    % stream power = "stream power" based bedrock erosion
    % both = 

% This code simulates fluvial erosion via two separate models, a
% transport-limited one with the Exner equation and a stream power law
% based detachment-limited one.

% CONSTANTS

rho = 1000; % water density, kg/m^3
rho_s = 2650; % sediment density, kg/m^3
g = 9.8; % acceleration due to gravity
D = 0.006; % sediment diameter
porosity = 0.35;

% SET-UP

L=100; %length in cells
dx=200; %cell length in meters
domain_length = L*dx;
x=(1:L);

z=.1*L-.1*x; %136*x.^(-0.3); % Elevations
z0=z; % Initial elevation pattern holder
deta=zeros(1,L); % Elevation change for each time-step

figure; hold;

%Preallocating
q_b_star=zeros(1,L);
q_b=zeros(1,L);
S=zeros(1,L);

dt=3600; % Time-step length, in seconds at bankfull discharge.

% Get discharge from Hack's law for drainage area as a function of
% downstream distance.

% A = ((L*dx-x)/C)^(5/7)

% From Rigon et al. (1996), Hack's law is:
% L = 1.4*A^0.6 with L in miles and A in square miles
% Converting to meters, L = 0.32*A^0.6
% Turning it around, A = (L/0.32)^(5/3)

A = (x*dx/0.32).^(5/3);

% Q = k A^n; log Q = log (k A^m); find slope
% (based on http://www.fgmorph.com/fg_7_6.php)
% F_0,x_0 = (50, 2) --metric--> (1.4, 5.2)
% F,x = (200, 10) --metric--> (5.7, 25.9)
% m = 0.86; k = 0.34
% With drainage area in square kilometers and discharge in m^3/s
% For both in meters, F_0,x_0 = (1.4, 5,200,000)
% m = 0.86, k = 2.35E-6

Q = 2.35E-6 * A.^.86;

% Use Leopold and Wolman's 0.4 exponent for downstream variation in
% bankfull depth with discharge. (Actually, right now, the coefficient and
% exponent are used based on the work of Wickert and Paola (in prep)
% because I (obviously) have it on-hand and L&W's 0.4 is approximate.

h = 0.46 * Q.^0.355;


% EVOLUTION

switch lower(type)
   
    
% Option 1: Exner equation
    case{'exner'}
      
for t=1:10000
    S(1) = (z(1)-z(2))/dx;
    S(L) = (z(end-1)-z(end))/dx;
    
    for x1 = 2:L-1
        S(x1) = (z(x1-1)-z(x1+1))/(2*dx);
    end
    
    for x1=1:L
        tau_b_star = rho*g*h(x1)*S(x1)/((rho_s-rho)*g*D);
        
        %MPM
        tau_c_star = 0.047;
        if tau_b_star > tau_c_star
            q_b_star(x1) = 8*(tau_b_star-tau_c_star).^1.5;
        else
            q_b_star(x1) = 0;
        end
        q_b(x1) = q_b_star(x1)*D*sqrt(((rho_s-rho)/rho)*g*D);
    end
    
    %Exner continuity
    
    for x1=2:L-1
        deta(x1) = (1/(1-porosity))*(q_b(x1-1)-q_b(x1+1))/(2*dx)*dt;
    end
    z=z+deta;
    
    if rem(t,100)==0
        plot(x*dx,z)
    end
    plot(x*dx,z0,'k')
end


% Option 2: Stream power law
    case{'stream power'}
        
for t=1:1000
    
    %Slope
    
    S(1) = (z(1)-z(2))/dx;
    S(L) = (z(end-1)-z(end))/dx;
    
    for x1 = 2:L-1
        S(x1) = (z(x1-1)-z(x1+1))/(2*dx);
    end
   
    %Stream power
    
    k = 0.00005;
    m = 1/3;
    n = 2/3;
        
    deta=-k*A.^m.*S.^n;
    deta(L)=0;
    z=z+deta;
    
end
        
    if rem(t,50)==0
        plot(x*dx,z)
    end
    plot(x*dx,z0,'k')
end