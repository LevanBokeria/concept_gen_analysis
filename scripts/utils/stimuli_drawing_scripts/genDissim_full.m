function [pairs, pair_indices, coords, nSesPairs] = genDissim_full

% This function creates the stimulus space, where the simuli are 16 points
% equally spaced from each other by 50 pixels, in 2 dimensions.

% Function will also create a matrix that specifies which pairs are
% compared for all 4 sessions of the task. session 3 and 4 are just repeats
% of session 1 and 2, but with mirrored location of pairs.

% Create the reward space
X = -75:50:75;
Y = X;

counter = 1;
for iX = 1:length(X)
    for iY = 1:length(Y)
        coords(counter,:) = [X(iX), Y(iY)];
        counter = counter + 1;
    end
end

% Calculate distances 
D = squareform(pdist(coords));

% Remove the ones that don't have the right distance.
% target_distances = [50,...
%                     100,...
%                     sqrt(50^2 + 50^2),...
%                     2*sqrt(50^2 + 50^2)];
% 
% 
% D(ismember(D(:),target_distances)) = 1;
% D(D ~= 1) = 0;

% Make the graph
G = graph(D);
pairs1 = G.Edges.EndNodes;

assert(size(pairs1,1) == 120,'Number of pairs selected is not right!\n');
% h = plot(G);
% h.XData = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
% h.YData = [repmat([1 2 3 4],1,4)];


%% Change coords to be in standard space, and shuffle them

nLoc = size(coords,1);

% Adjust coords by SF_dim scale
coords = coords + 175;

%% Create pairs for all sessions.
nPairs = size(pairs1,1);
nSesPairs = nPairs/2;

% Shuffle pairs along X
shuf_idx = randperm(nPairs);
shuf_idx = shuf_idx(1:nSesPairs);
pairs1(shuf_idx,:) = fliplr(pairs1(shuf_idx,:));

% Shuffle pairs along Y
pair_indices1 = (1:nPairs)';
pair_indices2 = (1:nPairs)';
pair_indices3 = (1:nPairs)';

shuffle_rows1 = randperm(nPairs);
shuffle_rows2 = randperm(nPairs);
shuffle_rows3 = randperm(nPairs);

pair_indices1 = pair_indices1(shuffle_rows1);
pair_indices2 = pair_indices2(shuffle_rows2);
pair_indices3 = pair_indices3(shuffle_rows3);

pairs1_shuffled = pairs1(shuffle_rows1,:);

%% Duplicate the pairs for 2nd round, in case subject has time
pairs2_shuffled = pairs1(shuffle_rows2,:);

% Switch the order of the 2nd round pairs
pairs2_shuffled = fliplr(pairs2_shuffled);

%% Make 3rd round
pairs3_shuffled = pairs1(shuffle_rows3,:);

% Concatenate the two matrices
pairs = [pairs1_shuffled; pairs2_shuffled; pairs3_shuffled];
pair_indices = [pair_indices1; pair_indices2; pair_indices3];

%% Double checking that pairs are shuffled correctly

% Check uniques
assert(isequal(unique(pairs2_shuffled(:)),unique(pairs1(:))),'Pairs messed up') ;
assert(isequal(sum(pairs2_shuffled(:)),sum(pairs1(:))),'Pairs messed up') ;

assert(isequal(pairs1(shuffle_rows1,:),pairs1_shuffled));
assert(isequal(pairs1(shuffle_rows2,:),fliplr(pairs2_shuffled)));

% Loop and check, I don't care!!!!
for iPair = 120:size(pairs,1)
    
    curr_pairs = pairs(iPair,:);
    
    match_pairs = pairs1(pair_indices(iPair),:);
    
    assert(isequal([1,1],ismember(curr_pairs,match_pairs)),'Pairs messed uuuup!');
    
end

end