%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numOriginalNodes = 24;
numAdditionalRenewableNodes = 8;     %this will be the number of new renewable energy sources we are adding to existing grid 
finalNumNodes = numOriginalNodes+numAdditionalRenewableNodes;   %total nodes when network is fully created
numEdges = 34;
numExtraEdges = 12;

%disp(A)
firstSum = sum(A);
avgDegree = sum(firstSum)/numOriginalNodes;                 %this is <k>
%/numNodes;
disp(sum(firstSum))
disp(avgDegree)

mySparseA = sparse(A);
%disp(mySparseA)

ChineseGraph = graph(A);
plot(ChineseGraph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nonrenewable portion of network successfuly generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start adding renewables (create different strategies for connections here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end adding renewables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start finding various properties of specific network (eigenvector
%centrality etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end finding various properties of specific network (eigenvector
%centrality etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start calculating stability using our stability metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We define the stability metric of our networks as being the time averaged
%difference between the normal frequency of 60 hertz and the average
%frequency of all oscillators at that time of the network 


%WHEN A RENEWABLE ENERGY SOURCE GOES OUT ITS FREQUENCY DROPS FROM 60 to 0
%IMMEDIATELY BECAUSE THERE IS NO INERTIA

%WHEN A NON RENEWABLE ENERGY SOURCE GOES OUT ITS FREQUENCY INCREASES THE
%LONGER IT IS DISCONNECTED (MUCH SMALLER PERTURBATION) WE WILL ASSUME THAT
%IT WILL BE FIXED ON AVERAGE WHEN ITS FREQUENCY HAS INCREASED FROM 60 TO
%60.5 HERTZ

probabilityRenewableOn = 0.3;
K = 3;      %coupling constant (this will only come into play when the nodes are connected on the grid)
dampingConstant = 0.7;
numTimeSteps = 3000;           %this is how long the simulations will run
changeTime = 0.01;
stabilityTimeVector = zeros(1, numTimeSteps);      %this will hold the stability metric of the network at each timestep (averaged at the end) 
singleOscillatorFrequency = zeros(1,numTimeSteps);
avgFrequencyFrom0 = zeros(1,numTimeSteps); 
coupledComponent1Timestep = zeros(1,numTimeSteps);
avgTheta = zeros(1,numTimeSteps);
stableWn = zeros(finalNumNodes,1);     %this is what I subtract from the resulting frequencies at each timestep when finding second norm

thetaVector = zeros(finalNumNodes,1);
dthetaVector = zeros(finalNumNodes,1);

WnVector = zeros(finalNumNodes,1);              %this is wN WHICH IS WHAT WE ARE CONCERNED ABOUT
dWnVector = zeros(finalNumNodes,1);

distFromTargetWn = zeros(finalNumNodes, 1);       %this will be used to keep track of how far each generator is from target frequency

powerGeneratedVector = randn(finalNumNodes,1);      %this is the power that each oscillator/generator consumes or produces
powerGeneratedVector = powerGeneratedVector-mean(powerGeneratedVector);      %this enforces that power in=power out on grid

dampingVector = dampingConstant*ones(finalNumNodes,1);

numTrials = 1;



for trial = 1:numTrials
    %must setup initial conditions of network (everything stable and good
    %initially at t=0
    
    %WnVector = stableWn;                   %PUT THIS ONE BACK JACKSON
    
    %WnVector = 0.01*randn(finalNumNodes,1);
    %WnVector(1)=50;
    WnVector = stableWn;
    
    omega = 0.01*randn(n,1);
    for i = 1:finalNumNodes          %here I setup the initial theta values of the generators (assuming uniform dist)
        initialAngle = rand()*2*pi;
        thetaVector(i) = initialAngle;       
    end
    
    
    iterator = 1;
    for t = 1:numTimeSteps
        
        
        dthetaVector = WnVector;            %change in theta = Wn
        dWnVector = powerGeneratedVector-dampingVector.*WnVector+ imag(exp(-1i*thetaVector).*(mySparseA*exp(1i*thetaVector)));      %change in frequency
        coupledComponent1Timestep = imag(exp(-1i*thetaVector).*(mySparseA*exp(1i*thetaVector)));
        
        %omegap = P - alpha.*omega + imag(exp(-1i*theta).*(a*exp(1i*theta)));
        
        
        
        
        
        
        
        %ABOVE I NEED TO CALCULATE dWnVector for this timestep
        
        WnVector = WnVector+changeTime*dWnVector;
        thetaVector = thetaVector+changeTime*dthetaVector;
        
        distFromTargetWn = WnVector-changeTime*stableWn;
        distFromTargetWn = distFromTargetWn.*powerGeneratedVector;  %here it is taken into consideration how much power a generator makes for how much it impacts stability
        stability = norm(distFromTargetWn-stableWn);
        stabilityTimeVector(t) = stability;
        
        singleOscillatorFrequency(t) = WnVector(1);
        
        avgFrequencyFrom0(t) = mean(WnVector.*abs(powerGeneratedVector));  %when the absolute value is there it shows
                                                                           %frequency of an oscillator*its power generated
        
          
    end
    
end

%plot(coupledComponent1Timestep)
%plot(stabilityTimeVector)
%plot(singleOscillatorFrequency)
plot(avgFrequencyFrom0)

%scatter(WnVector, powerGeneratedVector)

%plot stability against time
finalStabilityValue = mean(avgFrequencyFrom0);
finalStabilityValue = abs(finalStabilityValue);
disp("below is the final stability value for this network:")
disp(finalStabilityValue)




%STEPS TO DO:
%increase size of A to simulate adding in renewable energy sources

%run dynamics on the new matrix A where each row is an oscillator (kuramoto
%model)
%repeatedly alter A (aka cutting a connection in the grid and then test the
%stability)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





