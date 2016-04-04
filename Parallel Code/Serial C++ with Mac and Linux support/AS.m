% Even smaller AS data from Stanford SNAP
% David G. Anderson
% 2014
%

% read data
A  = textread('as19981229.txt');
A  = A(2:end, :);

% remove edges to self
if A(1,1)==A(1,2)
    A=A(2:end,:);
end
if A(end,1)==A(end,2)
    A=A(1:end-1,:);
end
for i=size(A,1)-1:-1:2
    if A(i,1)==A(i,2)
        A=[A(1:i-1,:);A(i+1:end,:)];
    end
end
A  = sort(A, 2);                        % double edges in this dataset %
nE = size(A, 1);
[d, I, AL] = unique([A(:,1); A(:,2)]);  % d(AL) = [A(:,1); A(:,2)] %
A2 = [AL(1:nE), AL(nE + 1: end)];
nN = max(max(A2));
AN = A2(:,1)*1e12 + A2(:,2);
[b1, b2] = sort(AN);
A2 = A2(b2, :);
A2 = A2(1:2:end, :);                    % remove duplicat rows %
nE = size(A2, 1);                       % update value %
S = sparse(A2(:,1), A2(:,2), ones(nE,1));
S(nN, nN) = 0;
S = S + S';
S = S + speye(nN);
% spy(S)

% form polar matrix
B = sparse(nE, nN);             % directed edge list
W = sparse(nE, nE);             % diagonal weight matrix
for iRow = 1:nE
    B(iRow, A2(iRow, 1)) = 1;
    B(iRow, A2(iRow, 2)) = -1;
    W(iRow, iRow) = 1;
end
Wh = 0 * W;
for iD = 1:nE
    Wh(iD, iD)  = sqrt(W(iD, iD));
end

% extract matrix U
% column selection in original problem equates 
% to row selection on U
[U, S, V] = svd(full(Wh * B));

% limit the problem to the data that is numerically nonzero
a = nnz(diag(S) > 1e-3);
U = U(:,1:a);

% parameters for column selection
k = nN - 1;                 % fixed number - not intended for asjustment
l = k + 1;                  % k < l <= nEs

% save for c++ implementation or use matlab implementation commented below
save('U.txt','U','-ascii','-double')

%{
% matlab column selectino algorithm (slow)
tic
pi = ucs(U(:, 1:a), a, 2 * a, false);
toc % 24 min
%}









