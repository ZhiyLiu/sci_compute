function f = randomRHS(dim)
    rand('seed', 23423);
    f = zeros(dim, dim);

%     d = repmat(dim, [1, dim]);
%     f(2:end-1, 2:end-1) = sin((x.*x +y.*y)./d);
    f(2:end-1, 2:end-1) = 1;%rand(dim-2, dim-2);
end