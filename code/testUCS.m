function testUCS
% Sparsification Examples for Unweighted Column Selection
% David G. Anderson
% 2015

% Two complete graphs connected by a single edge
fprintf('Example 1:\n');
fprintf('  Two complete graphs connected by a single edge\n')

% set coordinates
n = 9;
[coords, e] = CompleteConnected(n);
A = Adjacency(e);

% plot
gplot(A, coords); axis equal
set(findobj('Type','line'),'Color',[0.6 0.6 1]);
set(findobj('Type','line'),'LineWidth',0.5);
title('Connected Complete Graph')

% build matrix graph representation
[B, Wh] = GraphDecomp(e, 2*n);

% truncated SVD
[U, S, V] = svds(Wh * B, 2*n);

% unweighted column selection
k  = 2*n - 1;   % fixed number - not intended for asjustment
l  = 2 * k;     % k < l <= nEs
tic
idx = ucs(U(:, 1:end-1), k, l, false);
toc

% plot sparsifier
A2 = Adjacency(e(idx(1:k), :));
hold on
gplot(A2, coords, 'red');
title('Connected Complete Graph with Minimum Spanning Tree Sparsifier')
hold off
fprintf('Press Enter to continue\n');
pause

% dense sparsifier
fprintf('Example 2:\n');
fprintf('  Two complete graphs connected by a single edge with denser sparsifier\n')
gplot(A, coords); axis equal
set(findobj('Type','line'),'Color',[0.6 0.6 1]);
set(findobj('Type','line'),'LineWidth',0.5);
A2 = Adjacency(e(idx(1:l), :));
hold on
gplot(A2, coords, 'red');
title('Connected Complete Graph with Sparsifier with Oversampling')
hold off
fprintf('Press Enter to continue\n');
clear A A2
pause


% Random graph
fprintf('Example 3:\n');
fprintf('  Random graph\n')
n = 30;
[coords, e] = RandCoords(n, 10*n);
A = Adjacency(e);
gplot(A, coords); axis equal
set(findobj('Type','line'),'Color',[0.6 0.6 1]);
set(findobj('Type','line'),'LineWidth',0.5);
title('Random Graph')
[B, Wh] = GraphDecomp(e, n);
[U, S, V] = svds(Wh * B, n);
k = n-1;
l = 2*k;
tic
idx = ucs(U(:, 1:end-1), k, l, false);
toc
A2 = Adjacency(e(idx(1:k), :));
hold on
gplot(A2, coords, 'red');
title('Random Graph with Minimum Spanning Tree Sparsifier')
fprintf('Press Enter to continue\n');
pause

% denser graph
fprintf('Example 4:\n');
fprintf('  Random graph with denser sparsifier\n')
gplot(A, coords); axis equal
set(findobj('Type','line'),'Color',[0.6 0.6 1]);
set(findobj('Type','line'),'LineWidth',0.5);
A2 = Adjacency(e(idx(1:l), :));
hold on
gplot(A2, coords, 'red');
title('Random Graph with Sparsifier with Oversampling')
hold off
fprintf('End of examples\n');
clear A A2

end

function A = Adjacency(e)
% Returns adjacency matrix for a set of edges e
    n = length(unique([e(:,1); e(:,2)]));
    n = max(max(e));
    A = sparse(e(:,1), e(:,2), ones(size(e, 1), 1));
    A(n, n) = 0; % force matrix to be square
    A = A + A';
    A = A + speye(n);
end

function [B, Wh] = GraphDecomp(e, nN)
% Build signed edge-vertex incidence matrix
% and square root of weight matrix
    nE = size(e, 1);
    B = sparse(nE, nN);
    W = sparse(nE, nE);
    for iRow = 1:nE
        B(iRow, e(iRow, 1)) = 1;
        B(iRow, e(iRow, 2)) = -1;
        W(iRow, iRow) = 1;
    end
    Wh = 0 * W;
    for iD = 1:nE
        Wh(iD, iD)  = sqrt(W(iD, iD));
    end
end

function [coords, e] = CompleteConnected(n)
% Returns coordinates and edges for two complete graphs
% connected by a single edge
    coords=[];
    for i=0:n-1
        coords=[coords;[cos(i*pi*2/n)-1.5,sin(i*pi*2/n)]];
    end
    for i=0:n-1
        coords=[coords;[1.5-cos(i*pi*2/n),sin(i*pi*2/n)]];
    end
    % set edges
    e = nchoosek(1:n, 2);
    e = [e; nchoosek(1:n, 2)+n];
    e = [e; [1, n+1]];
end

function [coords, e] = RandCoords(n, l)
% Picks n random points in [0,1]x[0,1] and l random edges
    coords = rand(n, 2);
    e = zeros(l, 2);
    for i=1:l
        e(i, :) = randsample(n, 2, false);
        %e(i, :) = randsample(nchoosek(n,2), 1, false);
    end
end




