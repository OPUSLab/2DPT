function [J, h, cpy_idx, maximum_nbr] = sparse_J_h_input_copies(JJ,hb, W0, num_copies)

%% Be careful about processing the JJ matrix. Make it rounded to avoid asymmetry issues.

M = length(JJ);

%% Copy the initial matrix JJ to JJ2
JJ2 = JJ;

%% Calculate node degrees and number of copies required

node_degrees = max(degree(graph(JJ2)));  % Number of neighbors before sparsification

%num_copy_required = ceil(node_degrees / (max_neibr - 2)) % Copies including the original one

% if copies (including the original) goes beyond 2, then choose option2

option1_nbr = ceil(node_degrees / num_copies)+1; % Copies including the original one
option2_nbr = ceil(node_degrees / num_copies)+2; % Copies including the original one


if(num_copies == 2)
    num_bonds = option1_nbr-1;
else
    num_bonds = option2_nbr-2;
end

currentNode = M + 1;

%% Cell array to keep track of copies for each original node
copy_mapping = cell(M, 1);

for ii = 1:1:M
    for copy = 1:num_copies-1
        % Expand JJ2 by adding a column and row of zeros
        JJ2 = [JJ2 zeros(length(JJ2), 1)];
        JJ2 = [JJ2; zeros(1, length(JJ2))];

                % Keep track of the copy indices
        copy_mapping{ii} = [copy_mapping{ii}, currentNode];

        if(copy == 1)
            JJ2(ii, end) = W0;
            JJ2(end, ii) = W0;
        % elseif(copy == num_copy_required-1)
        %     JJ2(ii, end) = W0;
        %     JJ2(end, ii) = W0;
        %     JJ2(end-1, end) = W0;
        %     JJ2(end, end-1) = W0;
        else
            JJ2(end-1, end) = W0;
            JJ2(end, end-1) = W0;
        end

        % Find non-zero elements in the ii-th row
        [~, col, ~] = find(JJ2(ii, :));
        kkLim = min(length(col), num_bonds);

        % Redistribute connections
        for jj = 1:kkLim
            if jj <= length(col)
                JJ2(currentNode, col(jj)) = JJ2(ii, col(jj));
                JJ2(col(jj), currentNode) = JJ2(col(jj), ii);
                JJ2(ii, col(jj)) = 0;
                JJ2(col(jj), ii) = 0;
            end
        end
        currentNode = currentNode + 1;
    end
end

J=JJ2;

%% Finding the last index of every copies (including the original one)

all_neighbors = degree(graph(J))';
maximum_nbr= max(all_neighbors);

total_elements = length(all_neighbors);
last_indicies = (total_elements / num_copies) * (1:num_copies);

index_ends=last_indicies;

cpy_idx = cell2mat(copy_mapping);  % Copy indices as matrix

%% Distribute biases to original nodes and their copies

% Initialize the new bias vector for the expanded graph
hb_new = zeros(1, length(JJ2));

% Assign biases to original nodes and their copies
for ii = 1:M
    % Include the original node in the copies list
    all_copies = [ii, copy_mapping{ii}];
    num_copies = length(all_copies);
    
    % Distribute the bias equally
    hb_new(all_copies) = hb(ii) / num_copies;
end

h=hb_new;

end