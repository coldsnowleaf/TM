function [x, P] = FEKF(xtm1, u, P, zt, R, Q, AllRSS, Nap, K)
x = xtm1 + u; % Ex
P = P + R;  


xp = KNN(zt, AllRSS(1:Nap,:), AllRSS(Nap+1:Nap+2,:), K);

H = eye(2);
%H = [(zp - ztm1) ./ repmat(u(1),size(Q,1),1), (zp - ztm1) ./ repmat(u(2),size(Q,1),1)];% 20 * 2 用zt就没有错误
S = H*P*H' + Q; % 20*2 2*2 2*20 = 20*20 预测的方差

K = P*H'*pinv(S); %p*H'= 2*20 p*H'
x = x + K*(xp - H*x);
P = (eye(2) - K*H)*P;

end