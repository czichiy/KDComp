
%{

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

%}


%function T = CalculateT(KD,CvD, Nd)
function TD = CalculateTD(KDvec,CvD, Nd)


% component-wise form
% for i = 1 : Nd
%     for j = 1 : Nd
%          T(i,j) = ones(1,Nd)*(KD(:,:,(i-1)*Nd+j).*CvD')*ones(Nd,1); 
%     end
% end

% vectorized inner product (fast)
TD = reshape(KDvec'*reshape(CvD',[],1),Nd,Nd);

