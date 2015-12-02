% This is the main simulation file for grid search for. 
% 
% This code is for interest rate rule with INFLATION and UNEMPLOYMENT
%
% The logic of the code:
% 
% 1. Create intervals of parameters
% 2. Load intervals to dynare file and loop over values. Save dynare output from each iteration
% 3. Run welfare code with each of the output from dynare file
% 4. Find the highest welfare and corresponding parameter values


% 1. STEP

% define intervals for parameters in interest rate rule

inflBegin = 1.5; % 2.5
inflEnd   = 3; % 2.5,   3

uBegin    = 0.01;  % 0.01
uEnd      = 2;  % 1.5


%uBegin    = 0.01;  % 0.01
%uEnd      = 0.125;  % 1.5



StepSize  = 15;

% compute pairs

[ pairsParams ] = createIntervals( inflBegin, inflEnd, uBegin, uEnd, StepSize );

% save intervals for loading to dynare file

infl_Int = pairsParams(:,1);
u_Int    = pairsParams(:,2);

% delete previous versions if they exist

if ~isempty( ls('infl_Int.mat') )
    
    delete('infl_Int.mat')
    
end

if ~isempty( ls('u_Int.mat') )
    
    delete('u_Int.mat');
    
end


save infl_Int infl_Int
save u_Int u_Int

clear

% calculate welfare for each iteration

load paramValues

intervalLength = 1000;

storeWelfare = zeros(1,intervalLength);

for i = 1:intervalLength

	name              = sprintf('results_unemployment%d.mat',i); % set the output name

	inputFile         = load(name);                              % load the output

	storeWelfare(1,i) = welfareLoop( inputFile,paramValues );    % calculate and store the welfare

	i


end


% intervalLength = 1331;

% for i = 1:intervalLength

% 	name              = sprintf('results_unemployment%d.mat',i); % set the output name


% 	delete(name)


% end


% calculate welfare with different eta's

% load paramValues
% load eta_IntFig

% intervalLength = 1331;

% storeWelfareEta = zeros(1,intervalLength);

% for i = 1:intervalLength

% 	name              = sprintf('results_unemployment%d.mat',i); % set the output name

% 	inputFile         = load(name);                              % load the output

% 	storeWelfare(1,i) = welfareLoopEta( inputFile,paramValues,eta_IntFig(i,1) );    % calculate and store the welfare

% 	i


% end
