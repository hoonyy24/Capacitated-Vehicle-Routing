%% processHeuristicCVRPPlot_combined_deep_modified.m
% Improved ALNS-based CVRP solver complete code example (updateOperatorScores removed and route txt saving)
clc; clear; close all;
opengl('save','software');
set(groot, 'DefaultFigureVisible', 'off');

%% 1. Set up folders and instance settings
instanceFolder = 'C:\6002\instances';
solutionFolder = 'C:\6002\solutions';
outputFolder = 'C:\6002\Cluster';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

instanceFiles = dir(fullfile(instanceFolder, '*.vrp'));
numSelect = 1000; % Number of instances to process

Plot_data = struct('filename', {}, 'groupInfo', {}, 'instanceID', {}, 'totalCost', {}, ...
    'optimalCost', {}, 'costDifference', {}, 'solutionRoutes', {}, 'polarData', {}, ...
    'capacity', {});
nFiles = length(instanceFiles);
selectedIdx = randperm(nFiles, min(numSelect, nFiles));

%% 2. Processing loop for each instance
for i = 1:length(selectedIdx)
    k = selectedIdx(i);
    filePath = fullfile(instanceFolder, instanceFiles(k).name);
    fprintf('==== Processing instance file: %s ====\n', instanceFiles(k).name);
    
    % Read VRP instance (the readVRP function should be implemented separately)
    [coords, demand, depot, groupInfo, instanceID, cap] = readVRP(filePath);
    % The function returns 6 values in order:
    % cap = capacity
    demand = demand(:); % Ensure demand is a column vector
    if size(coords, 1) ~= length(demand)
        error('Inconsistent number of nodes: coords has %d rows, demand has %d elements in %s', ...
            size(coords, 1), length(demand), instanceFiles(k).name);
    end
    customers = [coords, demand];
    % Calculate number of vehicles required (excluding depot since its demand is 0)
    vehicleCount = ceil(sum(demand(2:end)) / cap);
    capacity = cap;
    
    % Pre-calculate the distance matrix (Euclidean distances between customer coordinates)
    distMat = squareform(pdist(customers(:,1:2)));
    fprintf('Group Info: %s, Instance ID: %s, Vehicle Count: %d, Capacity: %d\n', ...
        groupInfo, instanceID, vehicleCount, capacity);
    
    %% 3. Generate initial solution (based on Clarke-Wright Savings)
    initialSolution = generateInitialSolution_deep(customers, depot, vehicleCount, capacity);
    
    %% 4. ALNS Improvement Phase
    bestSolution = initialSolution;
    currentSolution = initialSolution;
    initializeOperatorScores();
    
    T = 1000;
    noImproveCount = 0;
    MaxIter = 100;
    DiversifyThreshold = 20;
    fraction = 0.2;  % Removal ratio
    
    for iter = 1:MaxIter
        [remOp, insOp] = selectOperators();
        [partialSol, removedList] = remOp(currentSolution, fraction, customers(:,1:2), customers(:,3), distMat);
        
        newSolution.routes = insOp(partialSol, removedList, customers(:,1:2), customers(:,3), capacity, distMat);
        newSolution.coords = customers(:,1:2);
        newSolution.cost = totalSolutionCost(newSolution.routes, distMat);
        
        % Apply light local search (2-opt) and inter-route swap
        newSolution = localSearchLight(newSolution);
        newSolution = interRouteSwap(newSolution, customers(:,3), capacity, distMat);
        
        if acceptSolution(newSolution, currentSolution, T)
            currentSolution = newSolution;
            if currentSolution.cost < bestSolution.cost
                bestSolution = currentSolution;
                noImproveCount = 0;
            else
                noImproveCount = noImproveCount + 1;
            end
        else
            noImproveCount = noImproveCount + 1;
        end
        
        T = coolDown(T, iter, MaxIter);
        
        % Diversification strategy
        if noImproveCount > DiversifyThreshold
            currentSolution = diversifySolution(bestSolution, customers, capacity);
            noImproveCount = 0;
        end
    end
    
    %% 5. Final improvement (apply VND/VNS)
    bestSolution = localSearchVND(bestSolution);
    totalCost = bestSolution.cost;
    
    %% 6. Compare with optimal solution (if available)
    solFileName = [instanceFiles(k).name(1:end-4) '.sol'];
    solutionFile = fullfile(solutionFolder, solFileName);
    if exist(solutionFile, 'file')
        [~, solCost, ~] = readSolution(solutionFile);  % readSolution function should be implemented separately
        fprintf('Optimal (solution file) total cost: %.2f\n', solCost);
        fprintf('Cost difference: %.2f\n', totalCost - solCost);
    else
        fprintf('No corresponding solution file found.\n');
        solCost = NaN;
    end
    
    % Save results in Plot_data structure
    Plot_data(i).filename = instanceFiles(k).name;
    Plot_data(i).groupInfo = groupInfo;
    Plot_data(i).instanceID = instanceID;
    Plot_data(i).capacity = capacity;
    Plot_data(i).totalCost = totalCost;
    Plot_data(i).optimalCost = solCost;
    Plot_data(i).costDifference = totalCost - solCost;
    Plot_data(i).solutionRoutes = bestSolution.routes;
    
    %% Save route information to a text file (Routes_filename.txt)
    validName = matlab.lang.makeValidName(instanceFiles(k).name);
    txtFileName = fullfile(outputFolder, ['Routes_' validName '.txt']);
    fid = fopen(txtFileName, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', txtFileName);
    end
    for r = 1:length(bestSolution.routes)
        route = bestSolution.routes{r};
        fprintf(fid, 'Route %d: ', r);
        fprintf(fid, '%d ', route);
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('Saved routes for %s to %s\n', instanceFiles(k).name, txtFileName);
    
    %% Calculate polar coordinate data
    numCustomers = size(customers,1) - 1;
    depotCoord = customers(1,1:2);
    polarData = zeros(numCustomers,2);
    for j = 1:numCustomers
        idx = j + 1;
        dx = customers(idx,1) - depotCoord(1);
        dy = customers(idx,2) - depotCoord(2);
        r = sqrt(dx^2 + dy^2);
        theta = atan2d(dy, dx);
        if theta < 0
            theta = theta + 360;
        end
        polarData(j,:) = [theta, log(r+1)];
    end
    Plot_data(i).polarData = polarData;
    
    %% Generate and save plots (Polar & Cartesian)
    fprintf('Generating Polar Plot...\n');
    figPolar = figure;
    polaraxes; hold on;
    colors = lines(numCustomers);
    for j = 1:numCustomers
        theta_rad = deg2rad(polarData(j,1));
        r = polarData(j,2);
        polarplot(theta_rad, r, 'o', 'Color', colors(j,:), 'MarkerFaceColor', colors(j,:));
    end
    title(['Polar Plot: ' instanceFiles(k).name ' (Group: ' groupInfo ', Inst: ' instanceID ')']);
    hold off;
    saveas(figPolar, fullfile(outputFolder, ['PolarPlot_' validName '.png']));
    close(figPolar);
    
    fprintf('Generating Cartesian Plot...\n');
    figCart = figure; hold on;
    for j = 1:length(bestSolution.routes)
        route = bestSolution.routes{j};
        if ~isempty(route)
            plot(customers(route,1), customers(route,2), '-o', 'LineWidth', 2);
        end
    end
    title(['Final Routes: ' instanceFiles(k).name ' (Group: ' groupInfo ', Inst: ' instanceID ')']);
    xlabel('X'); ylabel('Y'); grid on;
    hold off;
    saveas(figCart, fullfile(outputFolder, ['RoutePlot_' validName '.png']));
    close(figCart);
    
    fprintf('Saved plots for %s in folder %s\n', instanceFiles(k).name, outputFolder);
    fprintf('====================================\n');
end

%% Calculate Average Percentage Error
validIndices = ~isnan([Plot_data.optimalCost]);
if any(validIndices)
    totalCosts = [Plot_data(validIndices).totalCost];
    optimalCosts = [Plot_data(validIndices).optimalCost];
    errorRates = ((totalCosts - optimalCosts) ./ optimalCosts) * 100;
    meanErrorRate = mean(errorRates);
    fprintf('Average Percentage Error across all instances: %.2f%%\n', meanErrorRate);
else
    fprintf('No valid optimal costs found for error rate calculation.\n');
end

%% Subfunctions -------------------------------------------------------------------------
function sol = generateInitialSolution_deep(customers, capacity)
    % Generate initial solution using the Clarke-Wright Savings algorithm
    numNodes = size(customers,1);
    distMat = squareform(pdist(customers(:,1:2)));
    savings = [];
    % Node 1 is the depot, so start from node 2
    for i = 2:numNodes
        for j = i+1:numNodes
            saving = distMat(1,i) + distMat(1,j) - distMat(i,j); % Calculate saving
            savings = [savings; i, j, saving];
        end
    end
    savings = sortrows(savings, -3);
    % Initialize each customer as a separate route (depot -> customer -> depot)
    routes = cell(numNodes-1,1);
    for i = 2:numNodes
        routes{i-1} = [1, i, 1];
    end
    % Merge routes based on savings in descending order
    for s = 1:size(savings,1)
        i = savings(s,1);
        j = savings(s,2);
        route_i = find(cellfun(@(r) (length(r) >= 2) && (r(2)==i), routes));
        route_j = find(cellfun(@(r) (length(r) >= 2) && (r(2)==j), routes));
        if isempty(route_i) || isempty(route_j) || route_i == route_j
            continue;
        end
        newRoute = [routes{route_i}(1:end-1), routes{route_j}(2:end)];
        if routeDemand(newRoute, customers(:,3)) <= capacity
            routes{route_i} = newRoute;
            routes{route_j} = [];
        end
    end
    routes = routes(~cellfun(@isempty, routes));
    cost = totalSolutionCost(routes, distMat);
    
    sol.routes = routes;
    sol.cost = cost;
    sol.coords = customers(:,1:2);
end

function tot = routeDemand(route, demand)
    tot = sum(demand(route(2:end-1)));
end

function tot = totalSolutionCost(solution, distMat)
    tot = 0;
    for i = 1:length(solution)
        tot = tot + routeCost(solution{i}, distMat);
    end
end

function cost = routeCost(route, distMat)
    cost = 0;
    for i = 1:(length(route)-1)
        cost = cost + distMat(route(i), route(i+1));
    end
end

function initializeOperatorScores()
    fprintf('Operator scores initialized.\n');
end

function [remOp, insOp] = selectOperators()
    persistent remScores insScores remOps insOps;
    if isempty(remScores)
        remOps = {@randomRemoval, @demandRemoval, @distanceRemoval, @worstCostRemoval, @sequenceRemoval, @shawRemoval};
        insOps = {@greedyInsertion, @regretInsertion};
        remScores = ones(length(remOps),1);
        insScores = ones(length(insOps),1);
    end
    % Probability-based selection using Softmax
    remProb = exp(remScores) / sum(exp(remScores));
    insProb = exp(insScores) / sum(exp(insScores));
    
    remIdx = randsample(1:length(remOps), 1, true, remProb);
    insIdx = randsample(1:length(insOps), 1, true, insProb);
    
    remOp = remOps{remIdx};
    insOp = insOps{insIdx};
end

function Tnew = coolDown(T, iter, MaxIter)
    Tnew = T * (1 - iter/MaxIter);
    if Tnew < 1e-3
        Tnew = 1e-3;
    end
end

function newSol = diversifySolution(customers, capacity)
    fprintf('Diversifying solution...\n');
    % Example: Randomly rearrange all customers to create a new solution
    numCust = size(customers,1) - 1;
    customerIndices = randperm(numCust) + 1;  % Exclude depot
    newRoute = [1, customerIndices, 1];
    if routeDemand(newRoute, customers(:,3)) <= capacity
        routes = {newRoute};
    else
        routes = {};
        currentRoute = [1];
        currentLoad = 0;
        for i = 1:length(customerIndices)
            idx = customerIndices(i);
            if currentLoad + customers(idx,3) <= capacity
                currentRoute(end+1) = idx;
                currentLoad = currentLoad + customers(idx,3);
            else
                currentRoute(end+1) = 1;
                routes{end+1} = currentRoute;
                currentRoute = [1, idx];
                currentLoad = customers(idx,3);
            end
        end
        currentRoute(end+1) = 1;
        routes{end+1} = currentRoute;
    end
    distMat = squareform(pdist(customers(:,1:2)));
    newSol.routes = routes;
    newSol.cost = totalSolutionCost(routes, distMat);
    newSol.coords = customers(:,1:2);
end

function accept = acceptSolution(newSol, currentSol, T)
    deltaCost = newSol.cost - currentSol.cost;
    if deltaCost < 0
        accept = true;
    else
        prob = exp(-deltaCost / T);
        accept = (rand < prob);
    end
end

%% Removal Operators
function [partialSol, removedCustomers] = randomRemoval(solution, fraction, demand)
    removedCustomers = [];
    partialSol = solution.routes;
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        if length(rte) > 2
            custs = rte(2:end-1);
            numRemove = max(1, ceil(length(custs)*fraction));
            idxToRemove = weightedSampleNoReplacement(length(custs), numRemove, demand(custs));
            removed = custs(idxToRemove);
            removedCustomers = [removedCustomers, removed];
            rte(ismember(rte, removed)) = [];
            partialSol{i} = rte;
        end
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

function [partialSol, removedCustomers] = demandRemoval(solution, fraction, demand)
    allCust = [];
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        allCust = [allCust, rte(2:end-1)];
    end
    numRemove = max(1, ceil(length(allCust)*fraction));
    [~, sIdx] = sort(demand(allCust), 'descend');
    removedCustomers = allCust(sIdx(1:numRemove));
    partialSol = cell(size(solution.routes));
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        rte(ismember(rte, removedCustomers)) = [];
        partialSol{i} = rte;
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

function [partialSol, removedCustomers] = distanceRemoval(solution, fraction, distMat)
    depot = 1;
    allCust = [];
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        allCust = [allCust, rte(2:end-1)];
    end
    numRemove = max(1, ceil(length(allCust)*fraction));
    d2dep = distMat(depot, allCust);
    [~, sIdx] = sort(d2dep, 'descend');
    removedCustomers = allCust(sIdx(1:numRemove));
    partialSol = cell(size(solution.routes));
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        rte(ismember(rte, removedCustomers)) = [];
        partialSol{i} = rte;
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

function [partialSol, removedCustomers] = worstCostRemoval(solution, fraction, distMat)
    custInfo = [];
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        for j = 2:(length(rte)-1)
            c = rte(j);
            rcost = removalCost(rte, j, distMat);
            custInfo = [custInfo; c, rcost];
        end
    end
    numRemove = max(1, ceil(size(custInfo,1)*fraction));
    [~, idx] = sort(custInfo(:,2), 'descend');
    removedCustomers = custInfo(idx(1:numRemove),1)';
    partialSol = cell(size(solution.routes));
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        rte(ismember(rte, removedCustomers)) = [];
        partialSol{i} = rte;
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

function c = removalCost(route, pos, distMat)
    c = 0;
    if pos >= 2 && pos < length(route)
        before = route(pos-1);
        after = route(pos+1);
        original = distMat(before, route(pos)) + distMat(route(pos), after);
        direct   = distMat(before, after);
        c = original - direct;
    end
end

function [partialSol, removedCustomers] = sequenceRemoval(solution, fraction)
    removedCustomers = [];
    partialSol = cell(size(solution.routes));
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        if length(rte) > 3
            custs = rte(2:end-1);
            nCust = length(custs);
            numRemove = max(1, ceil(nCust * fraction));
            if nCust >= numRemove
                startIdx = randi(nCust - numRemove + 1);
                seq = custs(startIdx:startIdx+numRemove-1);
                removedCustomers = [removedCustomers, seq];
                rte(ismember(rte, seq)) = [];
            end
        end
        partialSol{i} = rte;
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

function [partialSol, removedCustomers] = shawRemoval(solution, fraction, distMat)
    allCust = [];
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        allCust = [allCust, rte(2:end-1)];
    end
    if isempty(allCust)
        removedCustomers = [];
        partialSol = solution.routes;
        return;
    end
    seed = allCust(randi(length(allCust)));
    similarities = zeros(1, length(allCust));
    for i = 1:length(allCust)
        c = allCust(i);
        similarities(i) = distMat(seed, c);
    end
    [~, sortedIdx] = sort(similarities, 'ascend');
    numRemove = max(1, ceil(length(allCust)*fraction));
    removedCustomers = allCust(sortedIdx(1:numRemove));
    partialSol = cell(size(solution.routes));
    for i = 1:length(solution.routes)
        rte = solution.routes{i};
        rte(ismember(rte, removedCustomers)) = [];
        partialSol{i} = rte;
    end
    partialSol = partialSol(~cellfun(@isempty, partialSol));
end

%% Insertion Operators
function newSolution = greedyInsertion(solution, removedCustomers, demand, capacity, distMat)
% removedCustomers: list of indices of removed customers
    newSolution = solution;
    for c = removedCustomers
        bestScore = Inf;
        bestIdx = -1;
        bestPos = -1;
        for i = 1:length(newSolution)
            rte = newSolution{i};
            if routeDemand(rte, demand) + demand(c) > capacity, continue; end
            for pos = 2:length(rte)
                newR = [rte(1:pos-1), c, rte(pos:end)];
                inc = routeCost(newR, distMat) - routeCost(rte, distMat);
                score = inc / (demand(c) + eps);
                if score < bestScore
                    bestScore = score;
                    bestIdx = i;
                    bestPos = pos;
                end
            end
        end
        if bestIdx == -1
            newSolution{end+1} = [1, c, 1];
        else
            rtmp = newSolution{bestIdx};
            newSolution{bestIdx} = [rtmp(1:bestPos-1), c, rtmp(bestPos:end)];
        end
    end
end

function newSolution = regretInsertion(solution, removedCustomers, demand, capacity, distMat)
    newSolution = solution;
    while ~isempty(removedCustomers)
        bestRegret = -Inf;
        bestCust = -1;
        bestDetail = [];
        for c = removedCustomers
            insCosts = [];
            details = [];
            for i = 1:length(newSolution)
                rte = newSolution{i};
                if routeDemand(rte, demand) + demand(c) > capacity, continue; end
                for pos = 2:length(rte)
                    newR = [rte(1:pos-1), c, rte(pos:end)];
                    inc = routeCost(newR, distMat) - routeCost(rte, distMat);
                    insCosts(end+1) = inc;
                    details = [details; i, pos, inc];
                end
            end
            if isempty(insCosts), continue; end
            insCosts = sort(insCosts);
            if length(insCosts) == 1
                regretVal = insCosts(1);
            else
                regretVal = insCosts(2) - insCosts(1);
            end
            if regretVal > bestRegret
                bestRegret = regretVal;
                bestCust = c;
                [~, minIdx] = min(details(:,3));
                bestDetail = details(minIdx,:);
            end
        end
        if bestCust == -1
            bestCust = removedCustomers(1);
            newSolution{end+1} = [1, bestCust, 1];
            removedCustomers(removedCustomers == bestCust) = [];
        else
            routeIdx = bestDetail(1);
            pos = bestDetail(2);
            rtmp = newSolution{routeIdx};
            newSolution{routeIdx} = [rtmp(1:pos-1), bestCust, rtmp(pos:end)];
            removedCustomers(removedCustomers == bestCust) = [];
        end
    end
end

%% Other Utilities
function idx = weightedSampleNoReplacement(n, k, weights)
    idx = zeros(k,1);
    avail = 1:n;
    if all(weights <= 0)
        for i = 1:k
            chosen = randi(length(avail));
            idx(i) = avail(chosen);
            avail(chosen) = [];
        end
        return;
    end
    for i = 1:k
        w = weights(avail);
        w = w + eps;
        w = w / sum(w);
        cdf = cumsum(w);
        r = rand;
        chosen_index = avail(find(cdf >= r, 1));
        if isempty(chosen_index)
            chosen_index = avail(end);
        end
        idx(i) = chosen_index;
        avail(avail == chosen_index) = [];
    end
end

function solClean = cleanRoutes(solution)
    solClean = {};
    for i = 1:length(solution)
        route = solution{i};
        if length(route) < 2, continue; end
        if length(route) == 2 && route(1) == 1 && route(2) == 1, continue; end
        newRoute = [route(1)];
        for j = 2:(length(route)-1)
            if route(j) ~= 1
                newRoute(end+1) = route(j);
            end
        end
        newRoute(end+1) = 1;
        if length(newRoute) < 3, continue; end
        solClean{end+1} = newRoute;
    end
end

function sol = localSearchLight(sol)
    distMat = squareform(pdist(sol.coords));
    for i = 1:length(sol.routes)
        rte = sol.routes{i};
        rte = twoOptRoute(rte, distMat);
        sol.routes{i} = rte;
    end
    sol.cost = totalSolutionCost(sol.routes, distMat);
end

function sol = localSearchVND(sol)
    distMat = squareform(pdist(sol.coords));
    for i = 1:length(sol.routes)
        rte = sol.routes{i};
        rte = vnsRoute(rte, distMat);
        sol.routes{i} = rte;
    end
    sol.cost = totalSolutionCost(sol.routes, distMat);
end

%% Local Search Subfunctions
function bestRoute = vnsRoute(route, distMat)
    bestRoute = route;
    bestCost = routeCost(bestRoute, distMat);
    neighborhoods = {@twoOptRoute, @threeOptRoute, @orOptRoute};
    k = 1;
    while k <= length(neighborhoods)
        newRoute = neighborhoods{k}(bestRoute, distMat);
        newCost = routeCost(newRoute, distMat);
        if newCost < bestCost
            bestRoute = newRoute;
            bestCost = newCost;
            k = 1;
        else
            k = k + 1;
        end
    end
end

function bestRoute = twoOptRoute(route, distMat)
    bestRoute = route;
    bestDist = routeCost(bestRoute, distMat);
    improved = true;
    while improved
        improved = false;
        for i = 2:(length(bestRoute)-2)
            for j = i+1:(length(bestRoute)-1)
                newR = [bestRoute(1:i-1), fliplr(bestRoute(i:j)), bestRoute(j+1:end)];
                newDist = routeCost(newR, distMat);
                if newDist < bestDist
                    bestRoute = newR;
                    bestDist = newDist;
                    improved = true;
                    break;
                end
            end
            if improved, break; end
        end
    end
end

function bestRoute = threeOptRoute(route, distMat)
    bestRoute = route;
    bestDist = routeCost(route, distMat);
    improved = true;
    while improved
        improved = false;
        for i = 2:(length(route)-4)
            for j = i+1:(length(route)-2)
                for k = j+1:(length(route)-1)
                    newR = threeOptSwap(route, i, j, k, distMat);
                    newDist = routeCost(newR, distMat);
                    if newDist < bestDist
                        bestRoute = newR;
                        bestDist = newDist;
                        improved = true;
                        break;
                    end
                end
                if improved, break; end
            end
            if improved, break; end
        end
    end
end

function newR = threeOptSwap(route, i, j, k, distMat)
    bestR = route;
    bestD = routeCost(route, distMat);
    cands = {};
    cands{1} = [route(1:i-1), fliplr(route(i:j-1)), fliplr(route(j:k-1)), route(k:end)];
    cands{2} = [route(1:i-1), fliplr(route(i:j-1)), route(j:k-1), route(k:end)];
    cands{3} = [route(1:i-1), route(i:j-1), fliplr(route(j:k-1)), route(k:end)];
    for n = 1:length(cands)
        cdist = routeCost(cands{n}, distMat);
        if cdist < bestD
            bestD = cdist;
            bestR = cands{n};
        end
    end
    newR = bestR;
end

function bestRoute = orOptRoute(route, distMat)
    bestRoute = route;
    bestDist = routeCost(route, distMat);
    improved = true;
    while improved
        improved = false;
        for segSize = 1:3
            for i = 2:(length(bestRoute)-segSize)
                seg = bestRoute(i:i+segSize-1);
                remain = [bestRoute(1:i-1), bestRoute(i+segSize:end)];
                for j = 2:length(remain)
                    newR = [remain(1:j-1), seg, remain(j:end)];
                    newD = routeCost(newR, distMat);
                    if newD < bestDist
                        bestRoute = newR;
                        bestDist = newD;
                        improved = true;
                        break;
                    end
                end
                if improved, break; end
            end
            if improved, break; end
        end
    end
end

function newSol = interRouteSwap(solution, demand, capacity, distMat)
    newSol = solution;
    improved = true;
    while improved
        improved = false;
        for i = 1:length(newSol.routes)
            for j = i+1:length(newSol.routes)
                route1 = newSol.routes{i};
                route2 = newSol.routes{j};
                for pos1 = 2:length(route1)-1
                    for pos2 = 2:length(route2)-1
                        cust1 = route1(pos1);
                        cust2 = route2(pos2);
                        newRoute1 = route1;
                        newRoute2 = route2;
                        newRoute1(pos1) = cust2;
                        newRoute2(pos2) = cust1;
                        if routeDemand(newRoute1, demand) <= capacity && routeDemand(newRoute2, demand) <= capacity
                            oldCost = routeCost(route1, distMat) + routeCost(route2, distMat);
                            newCost = routeCost(newRoute1, distMat) + routeCost(newRoute2, distMat);
                            if newCost < oldCost
                                newSol.routes{i} = newRoute1;
                                newSol.routes{j} = newRoute2;
                                improved = true;
                            end
                        end
                    end
                end
            end
        end
    end
    newSol.cost = totalSolutionCost(newSol.routes, distMat);
end
