[r,c] = size(Var);
uInner = Var(2:r-1,2:c-1);
% fInner = f(2:r-1, 2:c-1);
uStack = reshape(uInner', [],1);
% fStack= reshape(fInner', [],1);
fsuppose=An*uStack;
