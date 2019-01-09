function [U,T] = schurQRbasic(H)
% The input is an Hessenberg matrix
% Therefore, call the hessenberg reduction first and then call this routine
% Returns the schur factorization of the matrix A, such that,
% H = U*T*U'   OR    T = U'*H*U

maxit=500;
T=H;
U=eye(size(T)); % A is square matrix
for k = 1:maxit
    [Q,R] = qr(T);
    T = R*Q;
    U = U*Q;
end
end