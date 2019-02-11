function [x, P] = TEKF(xtm1, u, P, zt, R, Q, AllRSS, Nap, K)
x = xtm1 + u; % Ex
P = P + R;

%zp = net(x/StepSize);
%zp = zeros(size(zt));%meanRSSNorm;
%ztm1 = net(xtm1/StepSize);
%zp = zeros(size(zt));
%ztm1 = net((x-u)/StepSize);
%zp = zeros(size(zt))
zp = KNN(x, AllRSS(Nap+1:Nap+2,:),AllRSS(1:Nap,:), K);
ztm1 = KNN(xtm1, AllRSS(Nap+1:Nap+2,:),AllRSS(1:Nap,:), K);

% xt1 = KNN(ztm1, AllRSS(1:Nap, :), AllRSS(Nap+1:Nap+2, :), 3)
% Q = eye(2);
% K = P*pinv(P+Q);
% x = x + K*(xt - xp     
% P = (eye(2) - K)*P;


H = [(zp - ztm1) ./ repmat(u(1),size(Q,1),1), (zp - ztm1) ./ repmat(u(2),size(Q,1),1)];% 20 * 2 用zt就没有错误
S = H*P*H' + Q; % 20*2 2*2 2*20 = 20*20 预测的方差


K = P*H'*pinv(S); %p*H'= 2*20 p*H'
eu = K*(zt - zp) ./ u;
x = x + K*(zt - zp);
P = (eye(2) - K*H)*P;

end