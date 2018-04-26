function f = smoothRHS(dim)
    x = linspace(0,1, dim);
    y = linspace(0,1, dim);
    Sigma = [.25 .3; .3 1];
    [X,Y] = meshgrid(x,y);
    f = zeros(dim,dim);
    f = mvnpdf([X(:) Y(:)], 0.5, Sigma);
    f = reshape(f, dim, dim);
    f(1, :)  = 0;
    f(end, :) = 0;
    f(:, 1) = 0;
    f(:, end) = 0;
end