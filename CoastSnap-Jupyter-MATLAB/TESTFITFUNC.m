function [UV] = TESTFITFUNC(beta, xyz)

load('finalnonlinearfittest')

K = [fx 0 c0U;
    0 -fy c0V;
    0  0 1];

R = angles2R(beta(4), beta(5), beta(6));
IC = [eye(3) -beta(1:3)'];
P = K*R*IC;
P = P/P(3,4);

UV = P*[xyz'; ones(1,size(xyz,1))];
UV = UV./repmat(UV(3,:),3,1);

UV = [UV(1,:) UV(2,:)]';
