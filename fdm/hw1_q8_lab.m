w1 = [-1/12, -2/3, 0, 2/3, -1/12];
w2 = [-1/12, 4/3, -5/2, 4/3, -1/12];

% Wleft = [w1(1:2);w2(1:2)];
% Wright = [w1(3:4);w2(3:4)];
% 
% Uleft = [-1, -w1(4:5);-w2(3:5)];
% Uright = [1, -w1(1:2);-w2(1:3)];
% 
% exCoefLeft = Wleft\Uleft; % extra coefficients left
% exCoefRight = Wright\Uright; % extra coefficients right
% 
% fleft = Wleft\[0;fEval(xi)];
% fright = Wright\[0;fEval(xl)];
% 
% existingRowsLeft = -Uleft;
% existingRowsRight = -Uright;
% 
% newRowsLeft = existingRowsLeft + Wleft * 

h=1;N=3;
f = [1 1 1];

Aleft = zeros(3,6);
Aright = zeros(3,6);

Aleft(1,1:5) = w1;
Aleft(1,3) = -1;
Aleft(2,1:5) = w2;
Aleft(3,2:6) = w2;

Aright(1,2:6) = w1;
Aright(1,4) = 1;
Aright(2,2:6) = w2;
Aright(3,1:5) = w2;

bleft = [0 f(1) f(2)]'.*h^2;
bright = [0 f(N-1) f(N)]'.*h^2;

Bl = [Aleft,bleft];
Br = [Aright,bright];

Bl(2,:) = Bl(2,:) - Bl(1,:).*(Bl(2,1)/Bl(1,1));
Bl(3,:) = Bl(3,:) - Bl(2,:).*(Bl(3,2)/Bl(2,2));

Br(2,:) = Br(2,:) - Br(1,:).*(Br(2,6)/Br(1,6));
Br(3,:) = Br(3,:) - Br(2,:).*(Br(3,5)/Br(2,5));

