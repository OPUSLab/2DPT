function [connectivityMatrix] = read_Wishart_matrix(filename)

% Read the data from the file
data = readmatrix(filename);

% Determine the size of the connectivity matrix
% Assuming the indices are 1-based and the maximum index defines the size
maxIndex = max(data(:, 1:2), [], 'all');
connectivityMatrix = zeros(maxIndex, maxIndex);

% Fill the connectivity matrix with the values from the file
for i = 1:size(data, 1)
    row = data(i, 1)+1;
    col = data(i, 2)+1;
    value = data(i, 3);
    connectivityMatrix(row, col) = value;
    connectivityMatrix(col, row) = value;
end

% Display the connectivity matrix
disp('Connectivity Matrix:');
disp(connectivityMatrix);
end

