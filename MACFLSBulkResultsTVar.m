for boringstuff = 1:1       %not a for loop, just folds away comments
%A code which calls a linear multiple actuator cylinder code multiple
%times, to see how power varies when the row length varies.
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

%Any changes made to this code may also need made to
%MACFLSBulkResultsTConst.m, and vice versa!

%%%%%%This code is less polished compared to the others, being mostly copy-
%%%%%%pasted from the code I made for myself, so may constain errors or
%%%%%%inconsistencies, and will certainly contain notes to myself ^_^; 

end

%%

%Please note this code is less thoroughly checked others; it does
%seem to be working however :-)
for description_for_user = 1:1       %not a for loop, just folds away comments
%Code designed to be used to run MACFLSSolverVer.m multiple times, with a
%different number of turbines T each time.
%This paramater isn't varied in MACFLSBulkResultsTConst.m as changing T
%requires also changing the lengths of the passed TSR and tACW vectors.

%%%%%%MACFLSSolverVer.m depends on the following inputs, explained here:
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
%%% spacing         positive number                 spacing between adjacent turbines, measured in turbine diameters (=2R) (cosntant between all turbines) 
%%% Teval           vector length <=T containing
%%%                    natural numbers from 1-T     turbines from which results are extracted
%For avoidance of ambiguity, this is the format for variables being passed
%to MACFLSSolverVer.m to run it ONE time. If you want, for example, to use
%this code (MACFLSBulkResultsTVar) to consider a various different row
%lengths, then you will be setting T as a natural vector and passing each
%entry one at a time to the solver.

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
%(say) twelve T, then GeneratedPower will be a 12*length(Teval) matrix.

%%%%%%% ;
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
n = 16;         %(n must be divisible by 4)
U_0 = 10;
gamma = 10;
B = 2;
beta = 0;
R = 1;
c = 0.1;

%Set spacing between turbines.
%For multiple regularly-increasing spacings:
smin = 1.2;                                 %minimum spacing considered between turbines
smax = 4.2;                                 %maximum spacing considered between turbines
nS = 7;                                     %number of spacings to consider
spacing = smin:(smax-smin)/(nS-1):smax;     %array of all the considered spacings

%For specifically chosen spacings:
% spacing = [1.01:0.001:1.029 1.03:0.01:1.29 1.3:0.05:1.95 2:0.5:10];
% nS = length(spacing);

%If running the code only once:
% nS = 1; spacing = 3


%Set the number of turbines you want in each row that you're analysing
noTurbines = [2:4 6 12];
rT = length(noTurbines);            %number of rows you're considering
mT = max(noTurbines);               %number of turbines in longest row

% The following is a somewhat convoluted selection of options for various
% TSR and tACW combinations as the code cycles through varying T.
% The chosen method is to put the values into a matrix of size rT x mT,
% and then read the first noTurbines(i) entires for the i'th row.
% So eg, if noTurbines = [2:4 6 12], and the desired arrangement has the
% last turbine being of lower TSR than the others, then the resulting
% (rT=5 x mT=12) matrix is
% [4.8 4.5   0   0   0   0   0   0   0   0   0   0]     %two turbines
% [4.8 4.8 4.5   0   0   0   0   0   0   0   0   0]     %three turbines
% [4.8 4.8 4.8 4.5   0   0   0   0   0   0   0   0]     %four turbines
% [4.8 4.8 4.8 4.8 4.8 4.5   0   0   0   0   0   0]     %six turbines
% [4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.5]     %twelve turbines
% ..with the first 2,3,4,6,12 entries being used for the 1st, 2nd, 3rd,
% 4th, 5th row respectively.
% Both TSR and tACW  have an option for 'all the same', and an option for
% 'put arbritrary values in directly for each turbine/row combo'.
% Additional options have been provided for common other possibilities, eg
% alternating C/AC turbines or the above-mentioned "leading turbine has
% lower TSR than other turbines" (based on [ExpVal]'s oblique wind
% observation).
%%%%%%[ExpVal]
%
% Teval is kind of its own thing, so for now, I'm only allowing simpler
% options, namely analysis of the first and/or last and/or centre turbines
% in the rows, and analysis of every turbine in each row. Code exists to
% run evaluations for more general inputs, but you'll need to create the
% code to present such data yourself, potentially on a case-by-case basis.


%%%TSR:
%Same TSR for all turbines
TSR = 4.5*ones(rT,mT);

%End turbine has different TSR (put single entry at the front to make the
%leading turbine have the different TSR)
% TSR = zeros(rT,mT);
% for i = 1:rT
%     TSR(i,1:noTurbines(i)) = [4.8*ones(1,noTurbines(i)-1) 4.5];
% end

%For arbitrary TSRs (example based on noTurbines = [2:4 6 12])
% TSR = zeros(rT,mT);
% TSR(1,1:2) = [4.5,4.8];
% TSR(2,1:3) = [4.5,4.65,4.8];
% TSR(3,1:4) = [4.5,4.3,4.7,4.5];
% TSR(4,1:6) = [4.5,4.5,4.5,4.8,4.8,4.8];
% TSR(5,1:12) = [4.5 4.8*ones(1,11)];
TSR


%%%Rotation direction:
%All anticlockwise
% tACW = true(rT,mT);
%All clockwise
% tACW = false(rT,mT);                  

%First half one way, second half the other (can switch ceil/floor to alter 
%the centre turbine and true/false to change which half is AC or C)
% tACW = true(rT,mT);           %initialisation
% for i = 1:rT
%     tACW(i,1:noTurbines(i)) = [true(1,ceil(noTurbines(i)/2)) false(1,floor(noTurbines(i)/2))];
% end

%Alternating AC and C turbines (you can switch 'true' and 'false' for
%alternating C and AC turbines)
% tACW = true(rT,mT);
% for i = 1:rT
%     for j=2:2:noTurbines(i)     %works whether noTurbines(i) is even or odd
%         tACW(i,j) = false;
%     end
%     %%%tACW(2:2:noTurbines(i)) = false %???
% end

%For arbitrary rotation directions (example based on noTurbines = [2:4 6 12])
tACW = true(rT,mT);
tACW(1,1:2) = [false false];
tACW(2,1:3) = [false true false];
tACW(3,1:4) = [false true true false];
tACW(4,1:6) = [false true false false true false];
tACW(5,1:12) = [false true true true true false false true true true true false];

%another example
% tACW(1,1:6) = [true true true true true true];
% tACW(2,1:6) = [true true true false false false];
% tACW(3,1:6) = [false false false true true true];
% tACW(4,1:6) = [true false true false true false];
% tACW(5,1:6) = [true true false false true true];
% tACW(5,1:6) = [false false true true false false];

%For switching all turbines compared to last run
% tACW = not(tACW);


%%%Teval
%Setup
Teval = zeros(rT,mT);       %we filter out the nonzero entries later using nonzeros(Teval(iT,:))
allturb = false;            %used later
arbturb = false;            %used later

%Only the top/leftmost turbine each time
% Teval(:,1) = 1;
%Only the bottom/rightmost turbine each time
% Teval(:,1) = mT;

for i = 1:rT
    %Centre turbine    
    % Teval(i,1) = ceil(i/2);

    %A combination of two
    % Teval(i,1:2) = [1 ceil(noTurbines(i)/2)];

    %A combination of three
    % Teval(i,1:3) = [1 ceil(noTurbines(i)/2) noTurbines(i)];

    %All turbines in the row
    Teval(i,1:noTurbines(i)) = 1:noTurbines(i);
    allturb = true;

end                   

%Arbritrary input
%arbturb = true;
% Teval(1,1:2) = [1 2];
% Teval(2,1:1) = 1;
% Teval(3,1:3) = [1 3 4];
% Teval(4,1:2) = [1 6];
% Teval(5,1:12) = 1:12;

%Remember each row can't contain more entries than the max number of
%turbines you consider, shouldn't contain more entries than there are
%turbines in that run (unless you duplicate them for some reason), and
%cannot contain a number greater than the number of turbines in that run.
%So for the 2nd row of Teval with noTurbines = [2:4 6 12], Teval(2,1:3) = 
%[1 2 3] is fine; Teval(2,1:4) = [1 2 3 2] is valid but pointless;
%Teval(2,1:4) = [1 2 3 4]) is invalid; and Teval(2,1:15) =
%[1:3 1:3 1:3 1:3 1:3] is invalid (also pointless). Teval(2,1:3) = [2 3 1]
%seems to work according to brief testing.

%Now, we can pass nonzeros(Teval(iT,:))' to the solver (the transpose is
%needed to put the array in the right orientation).

%%%%%end user input required%%%%%



%%
P = zeros(nS,rT,mT);                        %Initialise power as function of spacing between turbines, which row you're considering (with a given # of turbines per row), and which turbine in each row you're analysing in that row (usually first, last, centre, or all of them; will have some zero entries if not the latter.)
Cp = zeros(nS,rT,mT);
%x_pertubation = zeros(nS,rT,n);             %initialise vector holding w_x, the normalised x-velocity pertubations
%y_pertubation = zeros(nS,rT,n);             %initialise vector holding w_y, the normalised y-velocity pertubations

for iT = 1:rT                               %for each row
    for iS = 1:nS                           %for each spacing
       %[soln,w_x,w_y,q,PowerGenerated,CpGenerated,u_tot,Re] = MACFLSSolverVer(n,U_0,gamma,T,B,beta,R,c,TSR,tACW,spac,Teval)         %for reference
        [~,w_x,w_y,q,PowerGenerated,CpGenerated,~,~] = MACFLSSolverVer(n,U_0,gamma,noTurbines(iT),B,beta,R,c,TSR(iT,1:noTurbines(iT)),tACW(iT,1:noTurbines(iT)),spacing(iS),nonzeros(Teval(iT,:))')          %Get (among other things) Power and Cp values at n theta-locations around the evaluated turbines
        %x_pertubation(kT,kS,:) = w_x(Teval,:);                                     %and only store results for the considered turbine.
        %y_pertubation(kT,kS,:) = w_y(Teval,:);                                     %and only store results for the considered turbine.
        P( iS,iT,1:length(nonzeros(Teval(iT,:))) ) = PowerGenerated(nonzeros(Teval(iT,:))');            %this might be transposed relative to how it should be
        %I think this indexing (ie PG(nonz(Teval)) ) works; moreover, I think it automatically handles duplicate evaluations, which is ideal. (Given Teval = [1 1 2], PG = [val1 val2] not [val1 val1 val2]; but, PG(Teval) = [PG(1) PG(1) PG(2)], so P(iT,iS,1:length(Teval)) is given length(Teval) entires, so we get no error.) 
        Cp( iS,iT,1:length(nonzeros(Teval(iT,:))) ) = CpGenerated(nonzeros(Teval(iT,:))');              %ditto above comments
    end
end

%%
%Now try to plot a 3D Power against spacing and row length/turbine position
%for the considered turbines
%This code doesn't work if Teval rows have leading zeros, but that error
%should never come up with the code as written.

if not(arbturb)             %don't do the below if Teval is arbitrary input; you can remove this
                            %if you want but it's probably worth just getting your plots manually
                            %and/or using the 'allturb' option after this if statement 
    
    %If only one turbine evaluated (commonly first, centre, or last turbine in the row)
    if isequal( Teval(:,1) , Teval(any(Teval,2),any(Teval,1)) )      %If matrix is empty except from a full first column, i.e. only one turbine evaluated in each row
        disp("if1")
        figure(4)
        %P = zeros(nS,rT,mT);                   %ordering reminder
        surf(noTurbines,spacing,P(:,:,1))       %a 1 here bc all the other columns of Teval are empty
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by first turbine at given spacing and row length')
        %title('Power generated by last turbine at given spacing and row length')
        %title('Power generated by centre turbine at given spacing and row length')
    
    
    %If two turbines evaluated (commonly first and last, or first/last and centre)
    elseif isequal( Teval(:,1:2) , Teval(any(Teval,2),any(Teval,1)) )       %If matrix is empty except from full first two columns, i.e. two turbines evaluated in each row
        disp("if2")
        %all of them on different plots
        figure(4)
        %P = zeros(nS,rT,mT);                   %ordering reminder
        surf(noTurbines,spacing,P(:,:,1),'FaceAlpha',0.5)       %a 1 here bc all the other columns of Teval are empty besides 2, done next
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by first/centre turbine at given spacing and row length')
        figure(5)
        surf(noTurbines,spacing,P(:,:,2),'FaceAlpha',0.5)       %a 2 here bc all the other columns of Teval are empty, done previous
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by centre/end turbine at given spacing and row length')
        
        %And/or all of them on same plot:
        figure(6)
        hold on
        %P = zeros(nS,rT,mT);                   %ordering reminder
        surf(noTurbines,spacing,P(:,:,1),'FaceAlpha',0.5)       %only 1 and 2 because all the other columns of Teval are empty
        % Phold = zeros(nS,rT);
        % for i = 1:rT
        %     Phold(:,i) = P(:,Teval(i,2));
        % end
        % surf(noTurbines,spacing,Phold,'FaceAlpha',0.5)           
        surf(noTurbines,spacing,P(:,:,2),'FaceAlpha',0.5)       %a 2 here bc all the other columns of Teval are empty, done previous
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by first/centre turbine at given spacing and row length')
    
    
    elseif isequal( Teval(:,1:3) , Teval(any(Teval,2),any(Teval,1)) )        %If matrix is empty except from full first three columns, i.e. three turbines evaluated in each row
        disp("if3")
        %all of them on different plots
        figure(4)
        %P = zeros(nS,rT,mT);                   %ordering reminder
        surf(noTurbines,spacing,P(:,:,1),'FaceAlpha',0.5)
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by first turbine at given spacing and row length')
        
        %title about for turbine centre or end (probably)
        figure(5)
        surf(noTurbines,spacing,P(:,:,2),'FaceAlpha',0.5)
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by centre turbine at given spacing and row length')
    
        %title about for turbine centre or end (probably)
        figure(6)
        surf(noTurbines,spacing,P(:,:,3),'FaceAlpha',0.5)
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by end turbine at given spacing and row length')
    
        %And/or all of them on same plot:
        %title about 'for turbine at 1+end or 1/end+centre'
        figure(7)
        hold on
        %P = zeros(nS,rT,mT);                   %ordering reminder
        surf(noTurbines,spacing,P(:,:,1),'FaceAlpha',0.5)
        surf(noTurbines,spacing,P(:,:,2),'FaceAlpha',0.5)
        surf(noTurbines,spacing,P(:,:,3),'FaceAlpha',0.5)
        
        % Phold2 = zeros(nS,rT);
        % for i = 1:rT
        %     Phold2(:,i) = P(:,Teval(i,2));
        % end
        % surf(noTurbines,spacing,Phold2,'FaceAlpha',0.5)           %only 1, 2, 3 because all the other columns of Teval are empty
        % 
        % Phold3 = zeros(nS,rT);
        % for i = 1:rT
        %     Phold3(:,i) = P(:,Teval(i,3));
        % end
        % surf(noTurbines,spacing,Phold2,'FaceAlpha',0.5)           %only 1, 2, 3 because all the other columns of Teval are empty
        % 
        ylabel('turbine spacing (D)')
        xlabel('number of turbines in row')
        zlabel('power')
        title('Power generated by first, centre, and last turbine at given spacing and row length')
    end

end
    
if allturb || arbturb                       %non-evaluated turbines will have P=0 (for arbturb)
% if allturb                                %can turn off arbturb if you want
    disp("if4")
    for i = 1:rT
        %for j = 1:noTurbines(i)
            figure(3+i) %: Teval(i,:)
            %P = zeros(nS,rT,mT);                   %ordering reminder
            surf(1:noTurbines(i),spacing,squeeze( P(:,i,1:noTurbines(i)) ),'FaceAlpha',0.5)       %Teval(i,:) = 1:noTurbines(i) so can use either 
            %This one makes more sense to plot spacing and the power of each turbine in the row, and have multiple figures for each row length
            %Contrast the previous where we plotted Power against spacing and row length, and used separate figures for which position we were looking at in the row (all the end turbines in the same figure, all the centre turbines in the same figure, etc)
            %Hence, P(:,i,:) rather than P(:,:,i)
            ylabel('turbine spacing (D)')
            xlabel('turbine position in row')
            zlabel('power')
            title('Power generated by each turbine in a row '+string(noTurbines(i))+' turbines, with given spacing between each turbine')
        %end
    end

end

%could also have something like 
% if arbturb
%     for i = 1:rT
%         figure(3+i)
%         surf(nonzeros(Teval(i,:))',spacing,squeeze( P(any(P,[2 3]),i,any(P,[1 2])) ),'FaceAlpha',0.5)
%     end
% end
%but it doesn't work as is, because squeeze(P(any(P,[2 3]),i,any(P,[1 2])))
%sometimes has zero arrays (maybe fix with rearrange(nonzeros(~),x,y) or do
%the P(any(P,x),any(P,y)) thing again?

%There are more things that can be done, but I want to get this code out
%onto github and then start looking at nonlinearities. Probably other
%general imporvements (eg, the case with arbturb = true) will gradually be
%added alongside the big improvements.

%Haven't checked but probably this is obsolete
% for i = 1:rT
%     figure(3+i)
%     surf(noTurbines,spacing,P(:,:,i))
%     ylabel('turbine spacing (D)')
%     xlabel('turbine count')
%     zlabel('power')
%     title('Power generated by turbine #'+string(Teval(i))+' at given spacing')
% 
%     % figure(4)
%     % hold on
%     % surf(noTurbines,spacing,P(:,:,1))
% end
% 
% 
% % %It's 4-dimensional time B-)
% % for i = 1:length(Teval)         %that ain't right
% %     figure(3+i)
% %     surf(noTurbines,spacing,P(:,:,1))       %that might not be right?
% %     xlabel('turbine spacing (D)')
% %     ylabel('turbine count')
% %     zlabel('power')
% %     title('Power generated by turbine #'+string(Teval(i))+' at given spacing')
% % end
% % 
% % for i = 1:length(Teval)
% %     figure(3+length(Teval)+i)
% %     surf(noTurbines,spacing,Cp(:,:,1))
% %     xlabel('turbine spacing (D)')
% %     ylabel('turbine count')
% %     zlabel('power')
% %     legend('Cp from turbine #'+string(Teval(i))+' at given spacing')
% % end
% % 
% % % figure(6)
% % % hold on
% % % %ldgXY=cell(1,length(Teval));
% % % plot(spacing,P(Teval,:))
% % % s = strings(length(Teval),1);
% % % for i = 1:length(Teval)
% % %     s(i) = 'Power generated by turbine #'+string(Teval(i))+' at given spacing';
% % % end
% % % legend(s)                %maybe a better way to do this, who cares though
% % % % 
% % % % for i = 1:T
% % % %     %ldgXY{i} = sprintf('average power of %d turbines',i);
% % % %     ldgXY{i} = sprintf('power of turbine %d in row of '+string(T)+' turbines',i);
% % % % end
% % % xlabel('turbine spacing (D)')
% % % ylabel('power')
% % % %legend(ldgXY,'Location','eastoutside')
% % 
% % figure(6)
% % hold on
% % %ldgXY=cell(1,rT);
% % for i = noTurbines
% %     plot(spacing,P(i,:))
% %     %ldgXY{i} = sprintf('average power of %d turbines',i);
% %     ldgXY{i} = sprintf('power of end turbine in row of %d turbines',i);
% % end
% % xlabel('turbine spacing (D)')
% % ylabel('power')
% % legend(ldgXY,'Location','eastoutside')
% % 
% % 
% % 
% % % figure(4)
% % % hold on;
% % % ldgXY=cell(1,nS);
% % % for k2 = 1:nS
% % %     plot(thetaval,x_pertubation(k2,:))                                      %plot the x-pertubations against theta on an x,y plot...
% % %     ldgXY{k2} = sprintf('spacing = %dD',spacing(k2));
% % % end
% % % legend(ldgXY,'Location','eastoutside')
% % 
% % % figure(5)
% % % ldgpolar=cell(1,nS);
% % % for k3 = 1:nS
% % %     polarplot(thetaval,abs(x_pertubation(k3,:)))                                 %...a polar plot...
% % %     hold on;                                                                    %(please note negatives are tricky with polar plots)
% % %     ldgpolar{k3} = sprintf('spacing = %dD',spacing(k3));
% % % end
% % % legend(ldgpolar)
% % 
% % % if nS>1
% % %     figure(6)
% % %     surf(thetaval,spacing,x_pertubation)                                        %...and a 3D plot.
% % %     xlabel('theta')
% % %     ylabel('turbine spacing (D)')
% % %     zlabel('w_x, x-velocity pertubation')
% % %
% % %     figure(7)
% % %     surf(thetaval,spacing,y_pertubation)                                        %...and a 3D plot.
% % %     xlabel('theta')
% % %     ylabel('turbine spacing (D)')
% % %     zlabel('w_y, y-velocity pertubation')
% % % end
% % 
% % % figure(8)
% % % plot(spacing,P)
