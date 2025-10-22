for boringstuff = 1:1       %not a for loop, just folds away comments
%A code which calls a linear multiple actuator cylinder code multiple
%times, where any of the variables can be altered except for the number of
%turbines in a row.
%An initial release of the planned full multiple actuator cylinder code.

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


%If you wish to receive updates when this code is updated; contact me with
%improvements, questions or suggestions regarding the code; or have found
%issues or errors, please feel free to contact me at
%angus.buck.2016@uni.strath.ac.uk (or, in the event my work email goes out
%of date, angus@alandish.co.uk).
%Future versions of this code may contain:
% - Handling of variables that change in tandem, eg U_0 and TSR pairs.
% - Showing variation of power as function of three variables via short
%   video output
%(Most improvements will occur in MACFLS.m and MACFLSSolverVer.m.)

%Changes made to this code may also need made to MACFLSBulkResultsTVar.m,
%and vice versa!

end

%%

for description_for_user = 1:1       %not a for loop, just folds away comments
%Code designed to be used to run MACFLSSolverVer.m multiple times, taking
%multiple n, U_0, gamma, B, beta, R, c, TSR layouts, tACW layouts, spacing,
%and Teval as input; and returning the corresponding
%[solution,w_x,w_y,q,PowerGenerated,CpGenerated,u_tot,Re] for each
%combination of variables. T must also be provided, but cannot be varied:
%if you wish to vary T, use MACFLSBulkResultsTVar.m instead.

%%%%%%MACFLSSolverVer.m depends on the following inputs, explained here:
%%% n               natural number                  number of evaluation points around each turbine (code is O(n^2); needs to be a multiple of 4)
%%% U_0             positive number                 ambient wind speed
%%% gamma           number                          incoming wind angle
%%% T               natural number                  number of turbines (code is O(T^2))
%%% B               natural number                  number of blades per turbine
%%% beta            number                          pitch angle of blades, in degrees
%%% R               natural number                  radius of turbines
%%% c               positive number                 aerofoil chord lengths
%%% TSR             positive vector, length T       tip-speed ratio of each turbine
%%% tACW            boolean vector, length T        rotation direction of each turbines (true = anticlockwise)
%%% spacing         positive number                 spacing between adjacent turbines, measured in turbine diameters (=2R) (cosntant between all turbines) 
%%% Teval           vector length <=T containing
%%%                    natural numbers from 1-T     turbines from which results are extracted
%For avoidance of ambiguity, this is the format for variables being passed
%to MACFLSSolverVer.m to run it ONE time. If you want, for example, to use
%this code (MACFLSBulkResultsTConst) to consider a turbine at twelve
%different U_0, then you will be setting U_0 as a positive vector and
%passing each entry one at a time to the solver.

%The following outputs are returned:
%%% w_x             matrix, size n*T                pertubation of the x-component of the wind velocity at each evaulated turbine point
%%% w_y             matrix, size n*T                    "       "      y-component  "      "       "       "       "       "       "
%%% q               matrix, size n*T                pressure jump/pressure exerted by turbine at each evaluated turbine point
%%% GeneratedPower  vector, length length(Teval)    power generated by each turbine
%%% CpGenerated     vector, length length(Teval)    coefficient of power for each turbine
%%% u_tot           matrix, size n*T                wind speed at each evaulated turbine point
%%% Re              matrix, size n*T                Reynolds number at each evaulated turbine point
%%% alpha           matrix, size n*T                angle of attack (deg) at each evaulated turbine point
%%% Cl              matrix, size n*T                coefficient of lift at each evaulated turbine point
%%% Cd              matrix, size n*T                coefficient of drag at each evaulated turbine point
%%% phi             matrix, size n*T                flow angle (rads) at each evaulated turbine point (equal to alpha*pi/180 if pitch = 0) 
% (note that swept area = 2R*H = 2R per unit height)
%Again, these are outputs from one run of MACFLSSolverVer.m. If you input
%(say) twelve U_0, then GeneratedPower will be a 12*length(Teval) matrix;
%if you additionally try three different pitches, GP will be a
%12*length(Teval)*3 multidimensional array, and so on.

%Relevant notes/protips:
% - squeeze(matrix) is used to remove dimensions of length 1 from a matrix.
%(For some reason, although selecting M(:,2) from a 2D matrix will return a
%vector, selecting M(:,2,5) from a 3D matrix will return a 3D matrix, where
%two of the dimensions are of length 1.)
% -  F5 will run the entire code. ctrl+enter will run the section of code
%currently selected, where a double percent/comment symbol ("%%") defines
%the start of a new section.
% - ctrl+R will comment out selected/highlighted lines of code, with a
%space between the % symbol and the text. shift+ctrl+R will instead remove
%a % symbol (and space if there is one).
%Based on this, commented lines without a space between the "%" symbol and
%the first word are just comments (like here); if there is a space, the
%comment is code that has been commented out, and is designed to be able
%to become code again with use of ctrl+shift+R.
end

%%

clearvars

%%%%%user input required%%%%%

%Set T: this is the only parameter which cannot be varied.
%If you want to get bulk results with rows of different lengths, please
%consult MACFLSBulkResultsTVariable.m.
T = 2;

%Set Teval
%In theory you could make this a variable too, but it would be entirely
%pointless, unless you're synchronising it with another variable (see
%later). 
Teval = 1:T;
% Teval = [2 3 5];


%Set n (all must be divisible by 4)
% n = 96;
% n = [16 32 64 128 256];
n = [8:8:32];
nn = length(n);

%Set U_0
U_0 = 8;
% U_0 = 8:0.5:12;
% U_0 = [6 8 9 9.5 10 10.5 11 12 15];
nU_0 = length(U_0);

%Set gamma
gamma = 0;
% gamma = 0:22.5:360;
% gamma = [-22.5 -10 0 10 22.5];
ngamma = length(gamma);
%gamma was added later than the others, so its implementation may be a bit
%less elegant

%Set B
B = 3;
% B = 2:3;
nB = length(B);

%Set beta
beta = 0;
% beta = -5 : 0.5 : 10;
% beta = [-2 -1 -0.5 0 0.5 1 2 3];
% beta = -1:1;
nbeta = length(beta);

%Set R
% R = 0.61;
R = 3:7;
% R = [1 10 100];
nR = length(R);

%Set c
c = 0.128;
% c = 0.4:0.2:1.6;
% c = [0.5 1 2];
nc = length(c);

%Set spacing between turbines.
% spacing = 0;
spacing = 1.5:0.5:5;
% spacing = [1.01:0.001:1.029 1.03:0.01:1.29 1.3:0.05:1.95 2:0.5:10];
% spacing = [1.05:0.05:1.5 1.6:0.1:2 2.25:0.25:10]
nS = length(spacing);


%TSR and tACW are a bit different, as they need an entire row to be passed,
%so varying them requires constructing a 2D array, not a 1D one.

%Set TSR
TSR = 2.7*ones(1,T)                 %no varying TSR layout between runs
%TSR= [2.7 2.7*ones(1,T-1)];        

% TSR(1,:) = 4.8*ones(1,T);
% TSR(2,:) = [4.5 4.8*ones(1,T-1)];
% TSR(3,:) = [4.8 5*ones(1,T-1)];            %varying TSR layout between runs

% TSR = [4.8*ones(1,T) ; 4.5 4.8*ones(1,T-1) ; 4.8 5*ones(1,T-1)]        %varying TSR layout between runs

nTSR = size(TSR,1);


%Set tACW
% tACW = false(1,T);                      %no varying TSR layout between runs
tACW = [true false];

% tACW(1,:) = false(1,T);                 %varying TSR layout between runs
% tACW(2,:) = true(1,T);
% tACW(3,:) = [false true false true false true];
% tACW(4,:) = [true false true false true false];
% tACW(5,:) = [true true true false false false];
% tACW(6,:) = not(tACW(5,:));

% tACW = [false(1,T); true(1,T) ; false true false true false true ; true false true false true false ; true true true false false false ; false false false true true true]

ntACW = size(tACW,1);


%we need one more user input later on, but this is fine for now

%%%%%end user input required%%%%%

%%

P = zeros(length(Teval),nn,nU_0,nB,nbeta,nR,nS,nc,nTSR,ntACW,ngamma);                               %initialise multidimensional array holding P, the Power output by each evaluated turbine for each run
Cp = zeros(length(Teval),nn,nU_0,nB,nbeta,nR,nS,nc,nTSR,ntACW,ngamma);                              %     "       "       "       "       "    Cp, the coefficient of power of each evaluated turbine for each run
% x_pertubation = zeros(T,nn,nU_0,nB,nbeta,nR,nS,nc,nTSR,ntACW,ngamma,max(n));                       %     "       "       "       "       "    w_x, the normalised x-velocity pertubations
% y_pertubation = zeros(T,nn,nU_0,nB,nbeta,nR,nS,nc,nTSR,ntACW,ngamma,max(n));                       %     "       "       "       "       "    w_y, the normalised y-velocity pertubations
%You can have length(Teval) or T as your dimension size in the first
%column; the former means you don't have rows of zeros for turbines you're
%not considering, the latter gives you an easier way to track which turbine
%in the row you're looking at. I've used T for the latter two as it makes
%extracting the data a bit easier in the for loop below ^_^;
%%%%%%This is probably not true (Teval(1:3) gives you T=1,6,12 easy
%enough); just work out if w_x(any(w_x,1),:) works.

%[solution,w_x,w_y,q,GeneratedPower,GeneratedCp,u_tot,Re,objectiveValue,exitflag,~] = MACFLSSolverVer(n,U_0,gamma,T,B,beta,R,c,TSR,tACW,spacing,Teval)

for iTSR = 1:nTSR
for itACW = 1:ntACW
    for in = 1:nn
    for iU_0 = 1:nU_0
    for igamma = 1:ngamma
    for iB = 1:nB
    for ibeta = 1:nbeta
    for iR = 1:nR
    for ic = 1:nc
    for iS = 1:nS
        %[soln,w_x,w_y,q,PowerGenerated,CpGenerated,u_tot,Re,objectiveValue,exitflag,~] = MACFLSSolverVer(n,U_0,gamma,T,B,beta,R,c,TSR,tACW,spacing,Teval)         %for reference
        
        [~,w_x,w_y,~,PowerGenerated,CpGenerated,~,~,objectiveValue,exitflag,~] = MACFLSSolverVer(n(in),U_0(iU_0),gamma(igamma),T,B(iB),beta(ibeta),R(iR),c(ic),TSR(iTSR,:),tACW(itACW,:),spacing(iS),Teval);
        
        % x_pertubation(:,in,iU_0,iB,ibeta,iR,iS,ic,iTSR,itACW,igamma,:) = w_x(Teval,:);
        % y_pertubation(:,in,iU_0,iB,ibeta,iR,iS,ic,iTSR,itACW,igamma,:) = w_y(Teval,:);
        %Note that for conditions where n(in) < max(n), the above rows will have zeros filling in entire from indices n(in)+1:max(n)
        P(:,in,iU_0,iB,ibeta,iR,iS,ic,iTSR,itACW,igamma) = PowerGenerated(any(PowerGenerated,1));
        Cp(:,in,iU_0,iB,ibeta,iR,iS,ic,iTSR,itACW,igamma) = CpGenerated(any(CpGenerated,1));
    end
    end
    end
    end
    end
    end
    end
    end
end
end

%This doesn't currently have the capability for having two variables alter
%at the same time. Eg, the 1.2kW Windspire validation has multiple wind
%speeds and *corresponding* TSRs (for a single turbine so TSR is a number,
%not a vector).
%If you wanted to do this, you could have:
 %Set corresponding pair
 % Pair(1,:) = [0.5:0.5:15.5]                     <<<--wind speed / corresponding TSRs vvv
 % Pair(2,:) = [0.01,0.57,0.91,1.87,1.15,1.36,1.34,1.60,2.00,2.16,2.12,2.09,2.12,2.18,2.18,2.21,2.22,2.20,2.18,2.15,2.10,2.04,1.95,1.84,1.70,1.59,1.39,1.02,0.86,0.78]
 % nPair = size(Pair,2)
%and relpace the U_0 and TSR for loops (or whichever other pair of for
%loops) with one for loop as below:
 % for iPair = 1:nPair
 %     MACFLSSolverVer(n(in),***Pair(1,iPair)***,gamma(igamma),T,B(iB),beta(ibeta),R(iR),c(ic),***Pair(2,Pair)***,tACW(itACW,:),spacing(iS),Teval);
 %     P = etc
 % end
%("***" only for emphasis)

%%

%This section will make plots of Power against each variable, for each
%turbine in Teval.

%This requires holding the other variables constant, which requires a
%'default' value of the other variables (if they are also varied).
%You can overwrite the below and choose your defaults individually;
%otherwise, the code below picks the middle index of each variable, thus
%taking the median value considered, which seems a reasonable way to
%select default values for the variables.
% dn     = ceil(nn/2);        dU_0   = ceil(nU_0/2);
% dB     = ceil(nB/2);        dbeta  = ceil(nbeta/2);
% dR     = ceil(nR/2);        dS     = ceil(nS/2);
% dc     = ceil(nc/2);        dTSR   = ceil(nTSR/2);
% dtACW  = ceil(ntACW/2);     dgamma = ceil(ngamma/2);

%Alternatively: This code gets the maximum power of all considered cases,
%and provides the indices where this occurs.
[M I] = max(P,[],"all");
P(I)
[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11] = ind2sub(size(P),I);
%ii = [i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11]
P(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11)

%Setting the default values for each variable to the one at which maximum
%power occurs, may provide a useful showcase of how each variable changes
%the power when everything else is kept optimal.
%Obviously, this is only for one turbine from Teval, and variables may be 
%correlated, so this isn't idealised information by any means, but the code
%is provided nonetheless.
dn     = i2;        dU_0   = i3;
dB     = i4;        dbeta  = i5;
dR     = i6;        dS     = i7;
dc     = i8;        dTSR   = i9;
dtACW  = i10;       dgamma = i11;


%%
%We plot {output} vs {variable} (say, Power vs U_0 or Cp vs TSR) for each
%evaluated turbine. These can be show either as a 3D plot or as stacked 2D
%plots.

% if size(P,7) ~= 1                       %if you're looking at the effect of different spacing (index 7)
%     % s = strings(1,length(Teval));
%     Pavg = zeros(1,nS);
%     figure(61)
%     hold on
%     for i = 1:nS
%         Pavg(i) = sum( squeeze(P(:,dn,dU_0,dB,dbeta,dR,i,dc,dTSR,dtACW,dgamma)))/2;
%     end
%     plot(spacing,Pavg)
%     xlabel('spacing between turbines, measured in turbine diameters')
%     ylabel('power')
%     % legend(s)
% end


%Stacked 2D plots section
for k = 1           %not a loop, I just want to be able to fold this code chunk away

if size(P,2) ~= 1                       %if you're looking at the effect of different n (index 2)
    s = strings(1,length(Teval));
    figure(22)
    hold on
    for i = 1:length(Teval)
        plot(n,squeeze(P(i,:,dU_0,dB,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('degree of code accuracy n')
    ylabel('power')
    legend(s)
end


if size(P,3) ~= 1                       %if you're looking at the effect of different U_0 (index 3)
    s = strings(1,length(Teval));
    figure(23)
    hold on
    for i = 1:length(Teval)
        plot(U_0,squeeze(P(i,dn,:,dB,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('ambient wind speed U_0')
    ylabel('power')
    legend(s)
end
 

if size(P,4) ~= 1                       %if you're looking at the effect of different B (index 4)
    s = strings(1,length(Teval));
    figure(24)
    hold on
    for i = 1:length(Teval)
        plot(B,squeeze(P(i,dn,dU_0,:,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('number of blades B')
    ylabel('power')
    legend(s)
end

if size(P,5) ~= 1                       %if you're looking at the effect of different beta (index 5)
    s = strings(1,length(Teval));
    figure(25)
    hold on
    for i = 1:length(Teval)
        plot(beta,squeeze(P(i,dn,dU_0,dB,:,dR,dS,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('pitch angle $beta$')
    ylabel('power')
    legend(s)
end

if size(P,6) ~= 1                       %if you're looking at the effect of different R (index 6)
    s = strings(1,length(Teval));
    figure(26)
    hold on
    for i = 1:length(Teval)
        plot(R,squeeze(P(i,dn,dU_0,dB,dbeta,:,dS,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('turbine radius R')
    ylabel('power')
    legend(s)
end

if size(P,7) ~= 1                       %if you're looking at the effect of different spacing (index 7)
    s = strings(1,length(Teval));
    figure(27)
    hold on
    for i = 1:length(Teval)
        plot(spacing,squeeze(P(i,dn,dU_0,dB,dbeta,dR,:,dc,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('spacing between turbines, measured in turbine diameters')
    ylabel('power')
    legend(s)
end

if size(P,8) ~= 1                       %if you're looking at the effect of different c (index 8)
    s = strings(1,length(Teval));
    figure(28)
    hold on
    for i = 1:length(Teval)
        plot(c,squeeze(P(i,dn,dU_0,dB,dbeta,dR,dS,:,dTSR,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('blade chord length')
    ylabel('power')
    legend(s)
end

if size(P,9) ~= 1                       %if you're looking at the effect of different TSR (index 9)
    s = strings(1,length(Teval));
    figure(29)
    hold on
    for i = 1:length(Teval)
        plot(1:size(TSR,1),squeeze(P(i,dn,dU_0,dB,dbeta,dR,dS,dc,:,dtACW,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('Tip speed ratio arrangement')
    ylabel('power')
    legend(s)
end

if size(P,10) ~= 1                       %if you're looking at the effect of different tACW (index 10)
    s = strings(1,length(Teval));
    figure(20)
    hold on
    for i = 1:length(Teval)
        plot(1:size(tACW,1),squeeze(P(i,dn,dU_0,dB,dbeta,dR,dS,dc,dTSR,:,dgamma)))
        s(i) = "Turbine #"+i;
    end
    xlabel('turbine rotation direction arrangement')
    ylabel('power')
    legend(s)
end

if size(P,11) ~= 1                       %if you're looking at the effect of different gamma (index 11)
    s = strings(1,length(Teval));
    figure(21)
    hold on
    for i = 1:length(Teval)
        plot(gamma,squeeze(P(i,dn,dU_0,dB,dbeta,dR,dS,dc,dTSR,dtACW,:)))
        s(i) = "Turbine #"+i;
    end
    xlabel('incoming wind angle')
    ylabel('power')
    legend(s)
end

end

%3D plots section
for k = 1           %not a loop, I just want to be able to fold this code chunk away
if length(Teval) ~= 1                       %can't have a 3D plot if you're only looking at one turbine
    
    if size(P,2) ~= 1                       %if you're looking at the effect of different n (index 2)   %if length(n)>1 would also work?
    figure(32)
    surf(n,Teval,squeeze(P(:,:,dU_0,dB,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
    xlabel('evaluation points (accuracy)')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,3) ~= 1                       %if you're looking at the effect of different U_0 (index 3)
    figure(33)
    surf(U_0,Teval,squeeze(P(:,dn,:,dB,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
    xlabel('incoming wind speed')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,4) ~= 1                       %if you're looking at the effect of different B (index 4)
    figure(34)
    surf(B,Teval,squeeze(P(:,dn,dU_0,:,dbeta,dR,dS,dc,dTSR,dtACW,dgamma)))
    xlabel('number of blades')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,5) ~= 1                       %if you're looking at the effect of different beta (index 5)
    figure(35)
    surf(beta,Teval,squeeze(P(:,dn,dU_0,dB,:,dR,dS,dc,dTSR,dtACW,dgamma)))
    xlabel('blade pitch angle')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,6) ~= 1                       %if you're looking at the effect of different R (index 6)
    figure(32)
    surf(R,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,:,dS,dc,dTSR,dtACW,dgamma)))
    xlabel('blade radius')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,7) ~= 1                       %if you're looking at the effect of different spacing (index 7)
    figure(37)
    surf(spacing,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,dR,:,dc,dTSR,dtACW,dgamma)))
    xlabel('spacing between turbines (measured in turbine diameters)')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,8) ~= 1                       %if you're looking at the effect of different c (index 8)
    figure(38)
    surf(c,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,dR,dS,:,dTSR,dtACW,dgamma)))
    xlabel('blade chord lengths')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,9) ~= 1                       %if you're looking at the effect of different TSR combinations (index 9)
    figure(39)
    surf(1:nTSR,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,dR,dS,dc,:,dtACW,dgamma)))
    xlabel('tip-speed ratio setup')
    ylabel('evaluated turbine')
    zlabel('power')
    end
    
    if size(P,10) ~= 1                       %if you're looking at the effect of different rotation direction configurations (index 10)
    figure(30)
    surf(1:ntACW,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,dR,dS,dc,dTSR,:,dgamma)))
    xlabel('rotation direction pattern')
    ylabel('evaluated turbine')
    zlabel('power')
    end

    if size(P,11) ~= 1                       %if you're looking at the effect of different gamma (index 11)
    figure(31)
    surf(gamma,Teval,squeeze(P(:,dn,dU_0,dB,dbeta,dR,dS,:,dTSR,dtACW,:)))
    xlabel('incoming wind angle')
    ylabel('evaluated turbine')
    zlabel('power')
    end

end
end


%Section for plotting two variables against each other within one turbine
%Example below is for varying n and spacing, but can be replaced with
%whichever variables you want to plot :-)
for i = 1:length(Teval)
    figure(40+i)            %figure(40+Teval(i)) might be helpful too
    surf(spacing,n,squeeze(P(i,:,dU_0,dB,dbeta,dR,:,dc,dTSR,dtACW)))
    %note surf(spacing,n,P) despite the order within P being n then spacing
    xlabel('separation (D)')
    ylabel('degree of accuracy')
    zlabel('Power (W)')
end

%There are of course other plots you can do! More may be included in future
%versions :-)

%Thank you for using my work and making it this far! :D
