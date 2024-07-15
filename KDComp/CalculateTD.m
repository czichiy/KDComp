% Example input
KD = rand(4, 4, 16); % Example 3D matrix
CvD = rand(4, 4);    % Example matrix
Nd = 4;              % Dimension size

% Reshape KD to a vector form suitable for the function
KDvec = reshape(KD, [], 1);

% Calculate TD
TD = CalculateTD(KDvec, CvD, Nd);

% Display the result
disp(TD);