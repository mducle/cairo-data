function [M, k, n, name, pname, limit] = bfof_structure(M0, x)
% Cairo pentagonal structure for Bi4Fe5O13F
%
% Syntax
%
% [M, k, n, name, pname, limit] = bfof_structure(M0, x)
%
% x: Input parameter [a0, a34, ac] of angles
%    a0: phase angle of 3-fold spins w.r.t. a-axis
%    a34: angle between the 3-fold and 4-fold spins
%    ac: angle between the 4-fold spins and the intermediate layer spins
%    angles are in radians
%
% M0: Magnitude of spins as a row vector
%
% M: Nx3 Spin component matrix
% k: propagation vector in rlu as a row vector (always [0.5 0.5 0.5])
% n: Normal vector of spin plane (always [0 0 1])
% name: name of function
% pname: parameter input name
% limit: upper/lower limits of parameters

if nargout < 4
    a0 = x(1); a34 = x(2); ac = x(3); b0 = x(4);
    S3 = [cos(a0); sin(a0); 0];
    S3p = [cos(a0+pi/2); sin(a0+pi/2); 0];
    T3 = [cos(b0); sin(b0); 0];
    T3p = [cos(b0+pi/2); sin(b0+pi/2); 0];
    S4 = [cos(a34); sin(a34); 0];
    S4p = [cos(a34+pi/2); sin(a34+pi/2); 0];
    Sc = [cos(ac); sin(ac); 0];
    Scp = [cos(ac+pi/2); sin(ac+pi/2); 0];
    M = [S4p S4 S4 S4p S4p S4 S4 S4p T3p -T3p -S3p S3p T3 -T3 S3 -S3 -Scp -Sc -Scp -Sc];
    M = [M -M -M M];
    k = [0 0 0];
    n = [0 0 1];
    if numel(M0) == 1
        M = M * M0;
    else
        M = M .* repmat(M0(:), 1, 3)';
    end
else
    name = 'Cairo magnetic structure for Bi4Fe5O13F';
    pname = {'Phi_0', 'Phi_34', 'Phi_c', 'Phi_1'};
    limit = [0 0 0 0;
             1 1 1 1] * (2 * pi);
    M = []; k = []; n = [];
end


