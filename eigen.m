function [V,D] = eigen(A)
% Returns the eigenvalues and eigenvectors of a general matrix A
% such that
% A = V'*D*V

[P,H] = hessenberg(A); % defined as A=PHP'
[Q,T] = schurQRbasic(H); % deefined as H=QTQ'

% Create the eigenvectors
Vtemp = P*Q;

% The eigenvalues are the diagonal elements of T
% Sort the eigenvalues and get the permutation vector
[Dv,I] = sort(diag(T));

% Extract the eigenvectors in the sorted order
V = Vtemp(:,I);
% Finally create the diagonal matrix for eigenvalues
D = diag(Dv);
end