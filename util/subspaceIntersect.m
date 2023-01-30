function X = subspaceIntersect(A,B,varargin)
%SUBSPACEINTERSECT Computes an orthonormal basis X for ran(A) intersected
%ran(B). SUBSPACEINTERSECT(A,B,tol) adds a tolerance threshold.
%   See Theorem 6.4.4 Golub, Matrix Computations

[theta, F, ~] = prinAngles(A,B);
if isempty(varargin)
    tol = 2*eps;
else
    tol = varargin{1};
end

angles = theta > (1 - tol);
X = F(:, angles);
end

