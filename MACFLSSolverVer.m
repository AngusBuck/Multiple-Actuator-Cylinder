%A linear multiple actuator cylinder code, an initial release of the
%planned full multiple actuator cylinder code.

%Produced by Angus Buck, PhD student at the University of Strathclyde,
%working under Dr Ed Hart and Prof. William Leithead (emeritus); formerly
%working under Peter Jamieson and Prof. J Michael Graham; and built
%following Andrew Ning (aning@byu.edu) in "Actuator cylinder theory for
%multiple vertical axis wind turbines" (2016) (doi:10.5194/wes-1-327-2016)
%in the Wind Enegy Science journal, itself derived from the original
%actuator cylinder model by Madsen (1982) (doi: 10.13140/RG.2.1.2512.3040)
%(link: dx.doi.org/10.13140/RG.2.1.2512.3040).

%COPYRIGHT
%The work in Ning's paper is provided under a Creative Commons 3.0 license
%(https://creativecommons.org/licenses/by/3.0).
%This code is derivative of, but not identical to, the work done by and
%software produced by Ning.
%Insofar as such activities do not infringe on Ning's CC 3.0 license,
%you are free to use, share, and edit this code, provided that credit is
%given to the author (Angus Buck, angus.buck.2016@uni.strath.ac.uk) and
%Andrew Ning.

%DISCLAIMER
%This code has not yet undergone rigorous testing and scrutiny
%to ensure it is without issue, and is unlikely to be considered without
%issue until further in my PhD, possibly not until I've defended my
%viva/undergone corrections, due starting Oct 2027.
%Additionally, even if the code is without error, it is currently only
%linearly accurate, so does not capture nonlinear effects, and in general
%is a low-to-medium-fidelity model. 
%Thus, the code is provided "as is", and no warranty is provided that it 
%provides accurate or meaningful results. Any user of the code accepts that
%the author has no responsibility for any losses or issues incurred as a
%result of errors or lack of accuracy this code contains, either in the
%software itself or the comments surrounding it; nor from any losses or
%issues resulting from use or abuse of the code if it is complete and
%functional.

%%%%%With the disclaimer made, a validation of the code is nonetheless
%provided alongside it, showing its accuracy compared to measured results
%and other's computaitonal results.


%If you wish to receive updates when this code is updated; contact me with
%improvements, questions or suggestions regarding the code; or have found
%issues or errors, please feel free to contact me at
%angus.buck.2016@uni.strath.ac.uk (or, in the event my work email goes out
%of date, angus@alandish.co.uk).
%Future versions of the code hope to contain:
% - Multiple rows and/or arbitrary turbine placement ('streamwise matrix of
%   turbines')
% - Nonlinear effects
% - Some level of 3D modelling, including stacked turbines ('vertical
%   matrix'), presumably multiple rows thereof (a tensor of matrices? :) )
%%%%%% - Tangential contribution to the power

%Any changes made to this code should (usually) also be made to MACFLS.m,
%and vice versa!


%%
function [solution,w_x,w_y,q,GeneratedPower,GeneratedCp,u_tot,Re,objectiveValue,exitflag,output,r,om] = MACFLSSolverVer(n,U_0,gamma,T,B,beta,R,c,TSR,tACW,spac,Teval)

%Eg: MACFLSSolverVer(36,10,4,3,-1,5,0.1,[4.5 4.8 4.8 4.8],[true false false true],1.5,[1 4])


%This code determines the pressure distribution around each turbine in a
%row of turbines; the resulting influence on the wind around each turbine;
%and the power generated per height by the turbine, as a function of
%several variables.
%Rather than being run itself and having variables defined by editing the
%code, it is designed to be called and have variables passed to it. In
%particular, this allows it to be called multiple times with multiple
%variables by the two MACFLBulkResults... files.

%The code depends on the following inputs, explained here:
%%% n               natural number                  number of evaluation points around each turbine (code is O(n^2); needs to be a multiple of four)
%%% U_0             positive number                 ambient wind speed
%%% gamma           number                          incoming wind angle
%%% T               natural number                  number of turbines (code is O(T^2))
%%% B               natural number                  number of blades per turbine
%%% beta            number                          pitch angle of blades, in degrees
%%% R               natural number                  radius of turbines
%%% c               positive number                 aerofoil chord lengths
%%% TSR             positive vector, length T       tip-speed ratio of each turbine
%%% tACW            boolean vector, length T        rotation direction of each turbines (true = anticlockwise)
%%% spac            positive number                 spacing between adjacent turbines, measured in turbine diameters (=2R) (cosntant between all turbines) 
%%% Teval           vector length <=T containing
%%%                    natural numbers from 1-T     turbines from which results are extracted

%The following outputs are returned:
%%% w_x             matrix, size n*T                pertubation of the x-component of the wind velocity at each evaulated turbine point
%%% w_y             matrix, size n*T                    "       "      y-component  "      "       "       "       "       "       "
%%% q               matrix, size n*T                pressure jump/pressure exerted by turbine at each evaluated turbine point
%%% GeneratedPower  vector, length length(Teval)    power generated by each turbine per unit height
%%% CpGenerated     vector, length length(Teval)    coefficient of power for each turbine
%%% u_tot           matrix, size n*T                wind speed at each evaulated turbine point
%%% Re              matrix, size n*T                Reynolds number at each evaulated turbine point
%(note that swept area = 2R*H = 2R per unit height)

%Relevant notes/protips:
%ctrl+R will comment out selected/highlighted lines of code, with a space
%between the % symbol and the text. shift+ctrl+R will instead remove a %
%symbol (and space if there is one).
%Based on this, commented lines without a space between the "%" symbol and
%the first word are just comments (like here); if there is a space, the
%comment is code that has been commented out, and is designed to be able
%to become code again with use of ctrl+shift+R.

%%

%Set turbine solidity
sol = c*B/R;
%There are different turbine solidity definitions, so be careful when
%comparing. Solidity here is mostly just used to tidy up later equations.

%Set the turbine centers' coordinates.
%In this code, the wind comes along the x-direction, and a turbine row
%that's perpendicular to the wind is along the y-axis.
%Note that the code mathematicallly always has wind coming from the left,
%so actually accounts for angled wind by altering the turbine layout.
%This means that gamma can also be used to model a staggered array of
%turbines.
turbineCentres = zeros(T,2);        %turbineCentres(tCI,:) = (x,y)-coordinates of centre of turbine tCI.
for tCI = 1:T
    turbineCentres(tCI,1) = ( spac*(tCI-1)*sin(gamma*pi/180) )*2*R;
    turbineCentres(tCI,2) = ( spac*(tCI-1)*cos(gamma*pi/180) )*2*R;
end

%Set the evaluation points, and angular distance between them.
deltatheta = 2*pi/n;
thetaval = (deltatheta/2 : deltatheta : 2*pi - deltatheta/2)' ;


%The velocity perturbation calculation can be separated into two
%components: the geometry, and the pressure distribution. This section
%calculates A_x and A_y, the geometry component.
%The geometry only depends on n, T, spacing, and gamma; so, it beings
%by checking if the calculation has been done before with these variables,
%and if it has, loads the results instead of going through all the
%integrations again. If they don't, the results are saved to a file once
%calculated, so the next run with the same variables takes much less time.
[x,y,D_x,Wake_x,A_y] = findDxWakexAy2DLinear(n,T,spac,gamma,turbineCentres/R,deltatheta,thetaval);
%Note that if you don't use this code as intended, it may give incorrect
%results. For example, if you have {3 turbines going from (0,0) to (6,8)}
%then switch to {3 turbines going from (3,4) to (-3,-4)} by changing the 
%code creating turbineCentres, then it will say the results are already
%saved and return those results despite the different geometry, because
%both are passing the variables gamma = 30, T=3, spac=5 (and the save file
%name doesn't read turbineCentres as this would require a huge complex
%filename).

%You can fairly easily modify turbineCentres to have arbitrary turbine
%locations, but you'd probably need to remove the variable-saving
%component (or modify it to save as "Layout1", "Layout2" etc).

%Note that if an error is found re. findDxWakexAy2DLinear.m, all the
%outputs saved by it may need deleted. Two earlier versions of the code had
%this issue, where the wrong 'turbineCentres' was passed, creating wrong
%Dx/Wakex/Ay outputs; then those wrongly-created files continued to be
%read once the code was changed, as the code was finding the outputs based
%on the passed T,n etc.

%Un-normalise (x,y) components
x = R*x;    y = R*y;
%Create A_x
A_x = D_x + Wake_x;

%%

%Get Cl, Cd, data for various alpha values.
%Data for a DU06W200 aerofoil from the github of Edgar Martinez-Ojeda1, an
%author on "Vertical-axis wind-turbine computations using a 2Dhybrid wake
%actuator-cylinder model" alongside Francisco Javier Solorio Ordaz1 and
%Mihir Sen (used for validation):
DU06W200CL = readmatrix("du06w200cl.csv");
DU06W200CD = readmatrix("du06w200cd.csv");
DU06W200AA = readmatrix("du06w200aa.csv");
DU06W200RE = readmatrix("du06w200re.csv");
alphavals = DU06W200AA;
ClData = DU06W200CL(:,6);       %%%%%%Column 1 had Cl = Cd????
CdData = DU06W200CD(:,6);
%"Fluid Dynamic Mechanisms..." (DOI:10.1016/j.renene.2016.08.015) also
%confirms that this is the aerofoil the Windspire uses.
%Plotting the data shows columns 5-6 as most similar to airfoiltools.com's
% lots (see below); presumably 5 is closest (as "0.1" (in reference)
%million = 100,000 (airfoiltools.com's value), but data for 6 is cleaner.

%Alternatively, data for a DU06W200 aerofoil @ Ncrit = 9, Re = 10^5, courtesy of http://airfoiltools.com.
%Not recommended, as alpha values often go higher than the range given by this site; the code won't crash if
%you do this, but the results will still be much less accurate. (I also don't really know what Ncrit is.)
%Code is nonetheless provided for using data provided in this format.
% DU06W200 = readmatrix("xf-du06-w-200-dt-50000");
% alphavals = DU06W200(:,1)';
% ClData = DU06W200(:,2)';
% CdData = DU06W200(:,3)';

%Alternatively, input data directly.
%Prone to typos, not recommended. This is NACA0018 data, left as an example.
% alphavals = -15.5:0.25:15.5;
% ClData = [-0.955,-1.0272,-1.0783,-1.1189,-1.1535,-1.1203,-1.1217,-1.1308,-1.1462,-1.1656,-1.1414,-1.1355,-1.1414,-1.1487,-1.121,-1.1195,-1.1138,-1.0953,-1.095,-1.0765,-1.0688,-1.0581,-1.0458,-1.0375,-1.0236,-1.0148,-1.0014,-0.9907,-0.979,-0.966,-0.9568,-0.941,-0.9288,-0.9162,-0.9009,-0.8871,-0.8738,-0.8569,-0.8409,-0.8252,-0.8096,-0.792,-0.7748,-0.7579,-0.7326,-0.6863,-0.6394,-0.5975,-0.5556,-0.5079,-0.4718,-0.4314,-0.386,-0.3503,-0.2985,-0.2642,-0.2136,-0.1834,-0.1421,-0.1085,-0.0739,-0.029,0,0.029,0.0739,0.1085,0.142,0.1834,0.2136,0.2641,0.2985,0.3503,0.386,0.4313,0.4718,0.5079,0.5555,0.5975,0.6394,0.6862,0.7326,0.7577,0.7746,0.7917,0.8094,0.8249,0.8406,0.8566,0.8736,0.8869,0.9007,0.9161,0.9286,0.9409,0.9567,0.9659,0.9789,0.9905,1.0013,1.0147,1.0236,1.0375,1.0458,1.0581,1.0689,1.0767,1.0951,1.0956,1.1142,1.1198,1.1215,1.1498,1.1419,1.1362,1.1426,1.1661,1.1468,1.1317,1.1232,1.1226,1.154,1.1193,1.0788,1.0273,0.9543];
% CdData = [0.08747,0.07547,0.06773,0.0621,0.05778,0.05664,0.05432,0.05169,0.04907,0.04675,0.04522,0.04388,0.04228,0.04066,0.03959,0.03838,0.03697,0.03611,0.03505,0.03402,0.03319,0.03218,0.03148,0.03057,0.02993,0.0291,0.02851,0.02777,0.02722,0.02658,0.02606,0.02552,0.02497,0.02454,0.0241,0.02367,0.02333,0.023,0.0227,0.02242,0.02221,0.02203,0.02186,0.02173,0.02166,0.02174,0.02182,0.02181,0.0218,0.02179,0.02171,0.02165,0.02156,0.02148,0.02133,0.02125,0.02102,0.02093,0.02074,0.02058,0.02059,0.02045,0.02047,0.02045,0.02059,0.02058,0.02074,0.02093,0.02102,0.02125,0.02133,0.02148,0.02156,0.02165,0.0217,0.02179,0.0218,0.02181,0.02182,0.02173,0.02165,0.02172,0.02186,0.02202,0.02221,0.02241,0.02269,0.023,0.02332,0.02366,0.0241,0.02454,0.02497,0.02551,0.02605,0.02658,0.02722,0.02777,0.02851,0.0291,0.02992,0.03056,0.03148,0.03218,0.03319,0.03402,0.03505,0.03611,0.03697,0.03838,0.0396,0.04067,0.04229,0.04388,0.04522,0.04677,0.04909,0.05171,0.05431,0.05659,0.05783,0.06216,0.06782,0.07564,0.08783];

%Alternatively, create polynomials and extract clean data from it.
%Not recommended, mostly legacy code for comparisons. Could feasibly be used to extrapolate data for greater alpha
%ranges, but in my testing, the range was so great that this was less accurate than just using alphamax or alphamin.
% ClFormula = polyfit(alphavals,ClData,3);                    %approximate ClData by a cubic formula
% CdFormula = polyfit(alphavals,CdData,2);                    %approximate CdData by a quadratic formula
% ClData = polyval(ClFormula,alphavals)
% CdData = polyval(CdFormula,alphavals)

%Note that these code lines work for the particular data files I've been
%using; you may need to edit this code and/or the "getCl()" and "getCd()"
%functions (see end of the file) to extract data from other sites.
%In particular, the later functions assume alphavals(1) is the lowest
%alpha value rather than the highest.

%%

%Initial guess for w = [w_x;w_y]
x0 = [-0.1*ones(n,T);0.05*ones(n/4,T);-0.05*ones(n/2,T); 0.05*ones(n/4,T)];         % w_x ~10% of U_0 and negative, w_y ~5% of U_0 in direction away from the VAWT centre.

%Combine A_x variable and A_y variable into one variable, so one solver can solve the equation
A = cat(1,A_x,A_y);

%a thing
% options.FunctionTolerance = 10^-6;

%Solve the nonlinear (transcendental?) system of equations.
% options = optimoptions("fsolve","Display","iter");          %I haven't investigated this thoroughly, might hold useful data
% [solution,objectiveValue] = fsolve(@(w) solveFlowPerturbation(A,U_0,w,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR),x0,options);       %Spits out a lot of text for a repeatedly-called function, so commented out.
[solution,objectiveValue,exitflag,output] = fsolve(@(w) solveFlowPerturbation(A,U_0,w,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR),x0);

function f = solveFlowPerturbation(A,U_0,w,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR)
    f = w - tensorprod(A,getq(U_0,w,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR),[2 4],[1 2]);
    %solves w = A*q(w), i.e. w-A*q(w) = 0, iteratively.

    %The called functions feed into each other as below:
    %                     /---> AoA --> Cl,Cd ----> q
    % U_n,U_t------> phi ------->------->------/
    %          \----->-----U_total^2----->----/
    %...where q is the turbines' pressure distribution.
end

if exitflag == -2
    om = output.message;
else
    om ="";
end

%Find w_x, w_y and q
w_x = solution(1:n,:)
w_y = solution(n+1:2*n,:)
q = getq(U_0,solution,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR)

%Find |u| and Re at each point around the turbine
u_tot = U_0*sqrt( (1+w_x).^2 + w_y.^2 )
Re = 1.225.*u_tot.*c./18.*10^6

%%
%A little section where I try to check/debug my fsolve function

r = sum(objectiveValue.^2,"all");     %This is the value that fsolve compares to sqrt(options.FunctionTolerance) = 10^-3

% figure(11)
% plot(sort(abs(objectiveValue)))
% ylabel("Aq - w")
% 
% figure(12)
% plot(log10(sort(abs(objectiveValue))))
% ylabel("log_{10}(Aq - w)")

%%
%A plot of wind velocity and components around each turbine
%This has been commented out, as it was slowing the code down and it's a
%bit weird having it when there's so many runs, but is kept for
%completeness

% figure(1)
% hold on
% quiver(x,y,ones(n,T),zeros(n,T),'off')          %ambient wind velocity
% quiver(x,y,w_x,w_y,'off')                       %induced wind velocity
% quiver(x,y,1+w_x,w_y,'off')                     %resultant wind velocity
% % quiver([-y -y -y],[x x x],[-w_y -w_y -zeros(n,T)],[1+w_x w_x ones(n,T)],'off')         %all three of the above
% %   %It would be great to scale these, but putting them all together makes
% %   %them all the same colour, and putting them apart means they scale rel-
% %   %-ative to themselves, not each other (so induced-velocity-arrows look
% %   %the same magitude as the other arrows bc they're scaled independantly)
% for tCI = 1:T
%     fplot(@(t) R*sin(t) + spac*(tCI-1)*sin(gamma*pi/180)*2*R, @(t) R*cos(t) + spac*(tCI-1)*cos(gamma*pi/180)*2*R)
%     text(turbineCentres(tCI,1),turbineCentres(tCI,2),'T'+string(tCI))
% end
% %Add title, labels, legend 
% title('Plot of each turbine and the wind influence')
% xlabel('velocity in x-direction')
% ylabel('velocity in y-direction')
% legend('initial wind direction, normalised by U_0','wind pertubation, normalised by U_0','resultant wind direction, normalised by U_0','Location','northwest')
% 
%If desired: rotate image 90 degrees (useful when |gamma| < 45deg as
%computer screens tend to be wider than they are tall).
% view(-(90),90)

%%

Cl = getCl(T,n,ClData,getclcdindex(U_0,solution,thetaval,alphavals,T,n,beta,tACW,TSR));
Cd = getCd(T,n,CdData,getclcdindex(U_0,solution,thetaval,alphavals,T,n,beta,tACW,TSR));
Phi = getPhi(U_0,solution,n,T,thetaval,tACW,TSR);
W2 = getW2(U_0,solution,n,T,thetaval,tACW,TSR);

PowerNormalisation = 0.5*1.225*U_0^3*2*R;     %To find Cp, power is normalised by 1/2 * rho * U_0^3 * swept area, ie the rectangle of 2R*(1) for per-height result.
%%%%%%^checkme

%Power information
WindPowerSeries = -(1.225*U_0^2)*U_0*( w_y.*cos(thetaval) - (1+w_x).*sin(thetaval) ) .*q *R*deltatheta;
WindPower = sum(WindPowerSeries,1)
WindCp = WindPower/PowerNormalisation
%Approximating the turbine as a continuous cylinder, all of which is extracting energy from the wind
%simultaneously, the power extracted from wind is rho*U_0^2*integral[u(theta)*q(theta)*R]d(theta).
%In our discretised model, this becomes sum{ rho*U_0^2*u(theta_i)*q(theta_i)*R*deltaTheta }.
%This is the upper limit of VAWT energy generation; unclear to me if this is the same as a Betz limit or not.

GeneratedPowerSeries = (TSR*U_0/R) *B/2/pi *R.*( Cl.*sin(Phi) - Cd.*cos(Phi) )*0.5*1.225.*W2*c *deltatheta;
%%%%%%This matches with [Ning 2016], and with [Madsen 1982] provided his Ft = ct*0.5*1.225.*W2*c
GeneratedPower = sum(GeneratedPowerSeries,1)
GeneratedCp = GeneratedPower/PowerNormalisation
%The power generated by a turbine blade at a given theta is omega*R*TangentialForce(theta) = TSR*U_0/R *R *[Cl*sin(phi) - Cd*cos(phi)] * 0.5*rho*W^2*c 
%This is then multiplied by the number of blades to determine the total power (energy generated per second) at any instant, and averaged by integrating over theta from 0->2pi and dividing by 2pi. 
%In our discretised model, this becomes omega*B/2pi*sum{R*TanForce*deltaTheta}.

%Plot power generated at each point around a revolution
%This has been commented out, as it was slowing the code down and it's a
%bit weird having it when there's so many runs, but is kept for
%completeness
% figure(3)
% hold on
% plot(thetaval,WindPowerSeries(:,Teval),thetaval,GeneratedPowerSeries(:,Teval))
% %%%%%%Recall that at each point, a single blade is technically generating GPS*2pi/deltaTheta/B energy per second. 
% %%%%%%(These plots show different things, but hopefully this graph + legend is still meaningful)
% title('Power generated at each point around a revolution')
% xlabel('theta')
% ylabel('Power generated')
% s = strings(1,2*length(Teval));
% for sI = 1:length(Teval)
%     s(sI) = 'Power extracted from the wind by turbine #'+string(Teval(sI))+' at given theta (all thetas generating power)';
%     s(length(Teval)+sI) = 'Power generated by turbine #'+string(Teval(sI))+' blade at given theta, normalised by total power generated in a cycle';
% end
% legend(s)

alpha = getAlpha(U_0,solution,n,T,thetaval,beta,tACW,TSR);


%%
%The various functions called, mostly as part of solveFlowPerturbation().
%The diagrams feed into each other as below:
%                     /---> AoA --> Cl,Cd ----> q (pressure distribution)
% U_n,V_t------> phi ------->------->------/
%          \----->-----U_total^2----->----/

function q = getq(U_0,w,thetaval,alphavals,T,n,ClData,CdData,sol,beta,tACW,TSR)
    q = sol/4/pi/(U_0^2) .* getCn(U_0,w,thetaval,alphavals,T,n,ClData,CdData,beta,tACW,TSR) .* getW2(U_0,w,n,T,thetaval,tACW,TSR);
end

function Cn = getCn(U_0,w,thetaval,alphavals,T,n,ClData,CdData,beta,tACW,TSR)
    Cn = getCl(T,n,ClData,getclcdindex(U_0,w,thetaval,alphavals,T,n,beta,tACW,TSR)).*cos(getPhi(U_0,w,n,T,thetaval,tACW,TSR)) + getCd(T,n,CdData,getclcdindex(U_0,w,thetaval,alphavals,T,n,beta,tACW,TSR)).*sin(getPhi(U_0,w,n,T,thetaval,tACW,TSR) );
end

function Vn = getVn(U_0,w,n,thetaval)
    Vn = U_0.*(1+w(1:n,:)).*sin(thetaval) - U_0.*w(n+1:2*n,:).*(cos(thetaval));
end

function Vt = getVt(U_0,w,n,T,thetaval,tACW,TSR)
    Vt = zeros(n,T);
    for I3 = 1:T
        if tACW(I3) == true
            Vt(:,I3) = U_0.*(1+w(1:n,I3)).*cos(thetaval) + U_0.*w(n+1:2*n,I3).*(sin(thetaval)) + TSR(I3)*U_0;        %=Omega(I)*R(I); I think this is always going to be equivalent as Omega is *defined* in terms of the other variables.
        else
            Vt(:,I3) = -U_0.*(1+w(1:n,I3)).*cos(thetaval) - U_0.*w(n+1:2*n,I3).*(sin(thetaval)) + TSR(I3)*U_0;       %=Omega(I)*R(I); I think this is always going to be equivalent as Omega is *defined* in terms of the other variables.
        end
    end
end

function W2 = getW2(U_0,w,n,T,thetaval,tACW,TSR)
    W2 = getVn(U_0,w,n,thetaval).^2 + getVt(U_0,w,n,T,thetaval,tACW,TSR).^2;
end

function phi = getPhi(U_0,w,n,T,thetaval,tACW,TSR)
    phi = atan(getVn(U_0,w,n,thetaval)./getVt(U_0,w,n,T,thetaval,tACW,TSR));
end

%assumes alpha(1) = minimum alpha and alpha(end) = maximum alpha
function clcdindex = getclcdindex(U_0,w,thetaval,alphavals,T,n,beta,tACW,TSR)
%possibly I can do this more efficiently, calling less often or putting it inside getCl/getCd
    clcdindex = zeros(n,T);
    % minValues(k1,k2) = zeros(n,T);
    % closestIndexes(k1,k2) = zeros(n,T);           %use this if you want to save the values I guess?
    holdalpha = getAlpha(U_0,w,n,T,thetaval,beta,tACW,TSR);
    for k1 = 1:n
        for k2 = 1:T
            if holdalpha(k1,k2) < alphavals(1) - 1
                disp('Warning: alpha more than 1 degree beyond lower limit.')
                clcdindex(k1,k2) = 1;
                %assumes alpha(1) = minimum alpha

            elseif holdalpha(k1,k2) > alphavals(end) + 1
                disp('Warning: alpha more than 1 degree beyond upper limit.')
                clcdindex(k1,k2) = length(alphavals);
                %assumes alpha(end) = maximum alpha
                
            else
                [~, closestIndex] = min(abs( alphavals - holdalpha(k1,k2) ));        %find the smallest difference between the calc'd alpha and the alphas available
                clcdindex(k1,k2) = closestIndex;
            end
        end
    end
end

function alpha = getAlpha(U_0,w,n,T,thetaval,beta,tACW,TSR)
    alpha = getPhi(U_0,w,n,T,thetaval,tACW,TSR)*180/pi - beta;     %we want alpha in degees, everything else in radians
end

%%%%%%This can probably be done without a loop
function Cl = getCl(T,n,ClData,clcdindex)
    Cl=zeros(n,T);
    for k1 = 1:n
        for k2 = 1:T
            Cl(k1,k2) = ClData(clcdindex(k1,k2));
        end
    end
end

function Cd = getCd(T,n,CdData,clcdindex)
    Cd=zeros(n,T);
    for k1 = 1:n
        for k2 = 1:T
            Cd(k1,k2) = CdData(clcdindex(k1,k2));
        end
    end
end

%This used to be a separate m-file; I can't remember if that had some
%benefit, but for now, it's back in the same m-file, but still a called
%function rather than intergrated into the code.
function [x,y,D_x,Wake_x,A_y] = findDxWakexAy2DLinear(n,T,spac,gamma,normalisedTurbineCentres,deltatheta,thetaval)

%This data is saved to a file with the following filename:
filename = "findDxWakexAxAyfiles\VariablesN"+string(n)+"T"+string(T)+"spac"+string(spac)+"gamma"+string(gamma)+".mat"
%Having "." in the filename should be fine. If you're worried, or know it 
%causes problems on your system, you could eg add a code line changing eg
%"1.5" to "1,5" or "1_5" etc?

if isfile(filename)
    disp("File with required variables found. Reading and returning required variables.")
    load(filename,"x","y","D_x","Wake_x","A_y")
else
    
    disp("File with required variables not found. Calculating and saving file with required variables.")
    
    %Initialise Variables
    x = zeros(n,T);                 %initialise absolute x-coordinates of evaluation points
    y = zeros(n,T);                 %     "     absolute y-coordinates      "       "   
    x_i = zeros(n,n,T,T);           %     "     relative x-coordinates      "       "   
    y_i = zeros(n,n,T,T);           %     "     relative y-coordinates      "       "   

    D_x = zeros(n,n,T,T);               %initialise matrix of influencing-on-influenced components (x-direction)
    Wake_x = zeros(n,n,T,T);            %     "       "       influencing's-wake-on-influenced components (x-direction)
    A_y = zeros(n,n,T,T);               %     "       "       influencing-on-influenced components (y-direction)
    %A_x                                %A_x = D_x + Wake_x doesn't need initialised
    
    %Calculate D_x,Wake_x,A_y
    for I = 1:T             %cycling through each turbine; index is over T
                            %'influenced turbine'

        for i = 1:n             %cycling through each position on each turbine; index is over n
                                %'influenced position on influened turbine'

            %The (x,y) values for each influenced point on each influenced turbine
            x(i,I) = -sin(thetaval(i)) + normalisedTurbineCentres(I,1);
            y(i,I) =  cos(thetaval(i)) + normalisedTurbineCentres(I,2);

            for J = 1:T         %cycling through each turbine's influence on the current turbine+position; index is over T
                                %'influencing turbine'

                for j = 1:n             %cycling through each [point on a given turbine]'s influence on the given point on the given turbine; index is over n
                                        %'influencing position on influencing turbine'  
                    
                    %Next, we find the relative positions of each influenced point compared to each influencing point
                    %(It's more intuitive not to have this, and just have x(i,I)-nTC(J) in the D_x (etc) formulae; but this is useful elsewhere and likely(?) more computationally efficient.)
                    x_i(i,j,I,J) = x(i,I) - normalisedTurbineCentres(J,1);
                    y_i(i,j,I,J) = y(i,I) - normalisedTurbineCentres(J,2);
                    %For example, for a turbine influencing itself, x_i(i,j,I,I) = x(i,I)-nTC(I,1) = -sin(thetaval(i))+{3}-{3} = -sin(thetaval(i)), like a turbine centred on (0,0).
                    %Or, for a turbine 5 normalised units away, x_i(i,j,I,J) = x(i,I)-nTC(J,1) = -sin(th(i))+{1}-{6} = -sin(th(i))-5, so the x-point is between 4 and 6 units away from the turbine centre 
                    

                    %Caclulate the coefficient matrix D_x
                    if I~=J
                        %if I=/=J, D_x(i,j,I,J) = full integral
                        D_x(i,j,I,J) = 1/(2*pi).*integral(@(theta) ((x_i(i,j,I,J) -10^-4 +sin(theta)).*sin(theta) - (y_i(i,j,I,J)-cos(theta)).*cos(theta))./((x_i(i,j,I,J) -10^-4 +sin(theta)).^2+(y_i(i,j,I,J)-cos(theta)).^2),thetaval(j)-deltatheta/2,thetaval(j)+deltatheta/2);
                        %I know I know, I hate this tiny 10^-4 offset too, but it stops the code giving orange integration errors
                        %I think it's right in that we're analysing on the just-upwind side of the turbine, but it will need checked near the top, particularly if we
                        %introduce wake expansion and the top/bottom of the circle is no longer the definition of the upwind/downwind evaluaiton curves; and of
                        %course, ideally there would be a more elegent way of handling this.
                        %I initially used 10^-5 because that was the smallest power-of-ten offset that worked, but I'm using 10^-4 now just in case 
                    else
                        %If I=J (i.e. finding the influence of the turbine on itself), we can calculate terms directly as below
                        %These matrices depend only on n, so it is possible to precompute and store if desired
                        if i~=j
                            D_x(i,j,I,J) = deltatheta/4/pi;
                        elseif i>n/2
                            D_x(i,j,I,J) = (1+1/n)/2;
                        else
                            D_x(i,j,I,J) = (-1+1/n)/2;
                        end
                    end

                    %Caclulate the coefficient matrix Wake_x

                    %If I=J, add wake terms for back half of turbine
                    %(We evaluate on the just-upwind sides of the disc, so the back is influenced by the front, but the front isn't influeced by the front and the back not by the back.)
                    %I believe this is fine due to the terms we chose for D_x with I=J,i=j.
                    if I==J
                        if i>n/2 && j==n+1-i
                            Wake_x(i,j,I,J) = -1;
                        end

                    %If I=/=J, add wake terms if I is behind J
                    elseif -1<=y_i(i,j,I,J) && y_i(i,j,I,J) <=1 && x_i(i,j,I,J) >=0 && x_i(i,j,I,J)^2 + y_i(i,j,I,J)^2 >=1

                        thetaJkTest = zeros(n/2,1);
                        for index = 1:n/2                                     %for the influencing turbine (loadform doesn't need to me symmetric for this to work)
                            thetaJkTest(index) = abs( acos(y_i(i,j,I,J)) - thetaval(index) );         %find the differences between {thetaJk from the influencing turbine} and the thetas corresponding to one of our existing panels (ngl I actually can't remember this for sure but it's something like that)
                                %note that acos(y_i) will be a +ve angle between 0 and pi, so
                                %this gives the index corresponding to the front of turbine J
                        end
                        thetaJkTest
                        [~,indexJk] = mink(thetaJkTest,1)                           %take the closest thetaval to thetaJk (we assume constant loading across panels so the value this angle produces is the same)
                        %w_x(I,i) = w_x(I,i) - p(indexJk) + p(n+1-indexJk);

                        if j==indexJk                       %for the component from the front of the influencing turbine
                            Wake_x(i,j,I,J) = -1;               % -q(acos(y))
                        elseif j==n-indexJk+1               %for the component from the back of the influencing turbine
                            Wake_x(i,j,I,J) = 1;                % +q(-acos(y))
                        end                                 %leading to the full -q(acos(y))+q(-acos(y)) term required for a point in the wake of J.

                    end

                    %Caclulate the coefficient matrix A_y
                    A_y(i,j,I,J) = 1/(2*pi).*integral(@(theta) ((x_i(i,j,I,J) -10^-4 +sin(theta)).*cos(theta) + (y_i(i,j,I,J)-cos(theta)).*sin(theta))./((x_i(i,j,I,J) -10^-4 +sin(theta)).^2+(y_i(i,j,I,J)-cos(theta)).^2),thetaval(j)-deltatheta/2,thetaval(j)+deltatheta/2);
                    
                end
            end
        end
    end
end

%save the variables for next time
save(filename,"x","y","D_x","Wake_x","A_y")

end

end