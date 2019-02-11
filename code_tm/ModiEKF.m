 function [x, P] = ModiEKF(xtm1, u, P, zt, ztm1, R, Q, Abtm1ReDB, AbtDB, Nap, K)
%StepSize = 20;
% R = eye(2)*10; 
% Q = eye(20); %WIFI
% StepSize = 20;
% al = 0.1;
% bl = 0.1;
% eSnag = 0.5 / StepSize;
% eSvar = 5;
% w0m = al/(1+al);
% wi = 1/(2*(1+al));
% w0c =  al/(1+al) + 1 - al^2 + bl;
% R = eye(2)*eSnag*eSnag;

%u = ptm1 - ptm2
x = xtm1 + u; % Ex
P = P + R;  
%u是期望测量

% x1 = xtm1 + u + al*repmat(eSnag,2,1);
% x2 = xtm1 + u - al*repmat(eSnag,2,1);
% Ex = w0m*x + wi*x1 + wi*x2;
% Sx = w0c*(x-Ex)*(x-Ex)' + wi*(x1-Ex)*(x1-Ex)' + wi*(x2-Ex)*(x2-Ex)';
%zp = zeros(size(zt));
% zp = KNN([(ztm1-RSSMean')./RSSVar'; repmat((u-NagMean')./NagVar', size(ztm1, 1)/2, 1) * ReScale], Abtm1ReDB', AbtDB', K).*RSSVar'+RSSMean';
% ztm1 = KNN([(ztm2-RSSMean')./RSSVar'; repmat((utm1-NagMean')./NagVar', size(ztm2, 1)/2, 1) * ReScale], Abtm1ReDB', AbtDB', K).*RSSVar'+RSSMean';
% zp = KNN([ztm1; u * ReScale], Abtm1ReDB', AbtDB', 1);
% ztm1 = KNN([ztm2; utm1 * ReScale], Abtm1ReDB', AbtDB', 1);
%zp = KNN([1/(RSSScaleMax-RSSScaleMin)*(ztm1-RSSScaleMin); (u - RePoMin)/(RePoMax-RePoMin)], Abtm1ReDB', AbtDB', K)*(RSSScaleMax-RSSScaleMin)+RSSScaleMin;
%ztm1 = KNN([1/(RSSScaleMax-RSSScaleMin)*(ztm2-RSSScaleMin); (utm1 - RePoMin)/(RePoMax-RePoMin)], Abtm1ReDB', AbtDB', K)*(RSSScaleMax-RSSScaleMin)+RSSScaleMin;
% zp = KNN([ztm1; u * ReScale], Abtm1ReDB', AbtDB', K);
% ztm1 = KNN([ztm2; utm1 * ReScale], Abtm1ReDB', AbtDB', K);
%  norm(zp-ztm1, 2)
%  norm(zp-zt, 2)
%  norm(ztm1-zt, 2)
 
zpind1 = KNNIndex([ztm1], Abtm1ReDB(:,1:Nap)',200);
zpind2 = KNNIndex([u], Abtm1ReDB(:,Nap+1:Nap+2)',200);
zpm1 = KNN([ztm1], Abtm1ReDB(:,1:Nap)', Abtm1ReDB(:,1:Nap)', K);
%zp = KNN([u], zpcad(1:2, :), zpcad(3:end, :), 5);
interAbtDB = AbtDB(intersect(zpind1, zpind2), :);
%zp = mean(interAbtDB(1:min(K,size(interAbtDB, 1)), :))';
zp = mean(interAbtDB)';

% zpm1cad = KNNRSS([ztm2], Abtm1ReDB(:,1:Nap)', [Abtm1ReDB(:,Nap+1:Nap+2), AbtDB]')
% zpm1 = KNNStep([utm1], zpm1cad(1:2, :), zpm1cad(3:end, :))

 
%zp = KNN([1/(RSSScaleMax-RSSScaleMin)*(ztm1-RSSScaleMin); (repmat(u, size(ztm1, 1)/2, 1)-RePoMin)/(RePoMax-RePoMin)], Abtm1ReDB', AbtDB', K)*(RSSScaleMax-RSSScaleMin)+RSSScaleMin;
%ztm1 = KNN([1/(RSSScaleMax-RSSScaleMin)*(ztm2-RSSScaleMin); (repmat(utm1, size(ztm1, 1)/2, 1)-RePoMin)/(RePoMax-RePoMin)], Abtm1ReDB', AbtDB', K)*(RSSScaleMax-RSSScaleMin)+RSSScaleMin;
% zp = net([ztm1', u']');
% ztm1 = net([ztm2', utm1']');




%zp = zeros(size(zt));
%zp = net([ztm1', (u/StepSize)']');%认为是期望
% ztm1 = net([ztm2', (utm1/StepSize)']');
% zp1 = net([ztm1', u1']');
% zp2 = net([ztm1', u2']');
% Ezp = w0m*zp + wi*zp1 + wi*zp2;
% Szp = w0c*(zp-Ezp)*(zp-Ezp)' + wi*(z1-Ezp)*(z1-Ezp)' + wi*(z2-Ezp)*(z2-Ezp)';
% Cxz = w0c*(x-Ex)*(zp-Ezp)' + wi*(x1-Ex)*(z1-Ezp)' + wi*(x2-Ex)*(z2-Ezp)';
% K = Cxz * pinv(Szp);
% x = x + K*(zp - zt);
% P = (eye(2) - K*H)*P;

% zt
% zp
% ztm1

H = [(zp - zpm1) ./ repmat(u(1),size(Q,1),1), (zp - zpm1) ./ repmat(u(2),size(Q,1),1)];% 20 * 2 用zt就没有错误
%H = [(ztm1 - ztm2) ./ repmat(xtm1(1)-xtm2(1),size(Q,1),1), (ztm1 - ztm2) ./ repmat(xtm1(2)-xtm2(2),size(Q,1),1)];% 20 * 2 用zt就没有错误
%这里的H类似C
S = H*P*H' + Q; % 20*2 2*2 2*20 = 20*20 预测的方差

K = P*H'*pinv(S); %p*H'= 2*20 p*H' 
eu = K*(zt - zp) ./ u;
x = x + K*(zt - zp);
P = (eye(2) - K*H)*P;


end