function [x, P] = MoveEKF(xtm1, u, P, zt, ztm1, R, Q, AbRSSDB, RePoDB, AbRSSMinusDB, Nap, K)
%StepSize = 20;
%x = xtm1 + u; % Ex
%P = P + R; 
x = xtm1 + u;
P = P + R;


%up = (KNN([(ztm1-RSSScaleMin)/(RSSScaleMax-RSSScaleMin); (zt-RSSScaleMin)/(RSSScaleMax-RSSScaleMin)], AbRSSDB', RePoDB', K))*(RePoMax-RePoMin)+RePoMin;
up = KNN([ztm1; zt], AbRSSDB', RePoDB', K);
%xp = KNN([zt-ztm1], AbRSSMinusDB', RePoDB', 3);
xp = xtm1 + up;


H = eye(2);
%H = [(zp - ztm1) ./ repmat(u(1),size(Q,1),1), (zp - ztm1) ./ repmat(u(2),size(Q,1),1)];% 20 * 2 用zt就没有错误
S = H*P*H' + Q; % 20*2 2*2 2*20 = 20*20 预测的方差

K = P*H'*pinv(S); %p*H'= 2*20 p*H'
x = x + K*(xp - H*x);
P = (eye(2) - K*H)*P;

end