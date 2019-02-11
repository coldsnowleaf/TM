clear
global Snag;
global ErrorMeanNg; 
global ErrorMeanEKF;
global ErrorMeanFP;
global ErrorMeanEKFHigh;
global ErrorMeanFPHigh;
global Svar;
global ErrorMeanArrayNg;
global ErrorMeanArrayEKF;
global ErrorMeanArrayFP;
global ErrorMeanArrayEKFHigh;
global ErrorMeanArrayFPHigh;
global Nap;
global Restart;
%global TestRSS;
global SnagTest;
global PathTest;
global ApLoc;
global PathTrain;
global NapPca;
global emp;
global ReScale;
global testrandp;

clear OnlineTest.mat
load OfflineTrain.mat
% load OfflineNNFP.mat
%load OfflineNN.mat
load OfflineNNMove.mat

%TestRSS = TrainRSS;
%Restart = 1;
% Snag = 0.5;
% Snag = 2;

PathNumTest = 100;
% if Restart == 1
%     for i = 1 : PathNumTest
%         TestRSS(i).Path(1,:) = rand(1,2) * Width;
%         TestRSS(i).Path(2,:) = rand(1,2) * Height;
%     end
% end
%Restart = 1;
if Restart == 1
    PathTest = [];
    for i = 1 : PathNumTest
        PathTest(i,:) = [rand(1,2) * Width, rand(1,2) * Height];
    end
end

for i = 1 : PathNumTest
    TestRSS(i).Path(1,:) = PathTest(i,1:2);
    TestRSS(i).Path(2,:) = PathTest(i,3:4);       
end

for i = 1 : PathNumTest
    %TestRSS(i).Path = TrainRSS(i).Path;
    temps = 0;
    for j = 2 : size(TestRSS(i).Path, 2)
        dis = norm(TestRSS(i).Path(:,j) - TestRSS(i).Path(:,j-1));
        samplePoints = round(dis/StepSize) + 1;
        for s = 1 : samplePoints
            temps = temps +1;
            offsetX = (TestRSS(i).Path(1, j) - TestRSS(i).Path(1, j-1)) / samplePoints * (s-1);
            offsetY = (TestRSS(i).Path(2, j) - TestRSS(i).Path(2, j-1)) / samplePoints * (s-1);

            for k = 1 : Nap
                d = pdist2(ApLoc(k,:), [TestRSS(i).Path(1, j-1) + offsetX, TestRSS(i).Path(2, j-1) + offsetY]);
                TestRSS(i).RSSs(k,temps) = St-S0-10*2*log(d/1)+normrnd(0,Svar);
                %TestRSS(i).RSSs(k,temps) = normrnd(0,Svar);
            end
            TestRSS(i).RSSs(k+1, temps) = TestRSS(i).Path(1, j-1) + offsetX;
            TestRSS(i).RSSs(k+2, temps) = TestRSS(i).Path(2, j-1) + offsetY;
        end
    end
    TestRSS(i).RSSsPca(1:NapPca, :) = (TestRSS(i).RSSs(1:Nap, :)'*mapping.M)';
end

% for i = 1 : PathNumTest
%     for j = 1 : fix(Nap*size(TestRSS(i).RSSs, 2)*emp)
%         w = fix(rand*size(TestRSS(i).RSSs, 2)) +1;
%         h = fix(rand*Nap)+1;
%         if w > size(TestRSS(i).RSSs, 2)
%             w = size(TestRSS(i).RSSs, 2);
%         end
%         if h > Nap
%             h = Nap;
%         end
%         TestRSS(i).RSSs(h,w) = -20;
%         TestRSS(i).RSSsPca(h,w) = -20;
%     end
% end

if Restart == 1
    testrandp = [];
    for i = 1 : PathNumTest
        testrandp(i).p = randperm(size(TestRSS(i).RSSs, 2)*Nap);
    end
end

for i = 1 : PathNumTest
    for j = 1 : fix(Nap*size(TestRSS(i).RSSs, 2)*emp)
%         w = fix(rand*size(TrainRSS(i).RSSs, 2)) +1;
%         h = fix(rand*Nap)+1;
        h = fix((testrandp(i).p(j)-1)/size(TestRSS(i).RSSs, 2))+1;
        %w = mod(testrandp(i).p(j),size(TestRSS(i).RSSs, 2))+1;
        w = testrandp(i).p(j)-(h-1)*size(TestRSS(i).RSSs, 2);
        if w > size(TestRSS(i).RSSs, 2)
            w = size(TestRSS(i).RSSs, 2);
        end
        if h > Nap
            h = Nap;
        end
        TestRSS(i).RSSs(h,w) = -20;
        TestRSS(i).RSSsPca(h,w) = -20;
    end
%     TrainRSS(i).RSSs(1:Nap, 1:size(TrainRSS(i).RSSs, 2)) = 0;
%     TrainRSS(i).RSSs(1,fix(size(TrainRSS(i).RSSs, 2)/2)+1) = 1000;
end

%Snag = 0.5;
% TestPointNum = 0;
% ErrorSumNg = 0;
% ErrorSumEKF = 0;
ErrorArrayNg =[];
ErrorArrayEKF = [];


if Restart == 1
    for i = 1 : PathNumTest
        for j = 1 : size(TestRSS(i).RSSs, 2) 
            SnagTest(i,j).Norm = [normrnd(0, Snag * StepSize), normrnd(0, Snag * StepSize)];
        end
    end
end
% PathNg = [];
% PathEKF = [];


% for K = 10: 10
    K = 3;
% for rq = -3: 1: 3%1: 2 %RQ
    rq = 0;
    %rq = 0;
    ErrorSumNg = 0;
    ErrorSumEKF = 0;
    ErrorSumFP = 0;
    ErrorSumEKFHigh = 0;
    ErrorSumFPHigh = 0;
    TestPointNum = 0;
    ErrorArrayNg = [];
    ErrorArrayEKF = [];
    ErrorArrayFP = [];
    ErrorArrayEKFHigh = [];
    ErrorArrayFPHigh = [];
    eu1a = [];
    eu2a = [];
    
    for i = 1 : PathNumTest
        PoNg = [];
        PoEKF = [];
        PoTrue = [];
        PoFP = [];
        PoEKFHigh = [];
        PoFPHigh = [];
        P = eye(2);
        PHigh = eye(2);
        PFP = eye(2);
        PFPHigh = eye(2);
        PoTrue = TestRSS(i).RSSs(Nap+1:Nap+2,:);
        if size(PoTrue, 2) >= 2
            PoNg(:,1:2) = PoTrue(:,1:2);
            %PoNg(:,2) = PoTrue(:,2);
            PoEKF(:,1) = PoNg(:,1); 
            PoEKFHigh(:,1) = PoNg(:,1);
            PoFP(:,1) =  PoNg(:,1);
            PoFPHigh(:,1) =  PoNg(:,1);
            %oneStepMeasureBefo = (PoTrue(:,2)-PoTrue(:,1)) + SnagTest(i,2).Norm';
            for j = 2: size(PoTrue,2) 
                oneStepMeasure = (PoTrue(:,j)-PoTrue(:,j-1)) + SnagTest(i,j).Norm';
                %oneStepMeasure = (PoTrue(:,j)-PoTrue(:,j-1)) + [normrnd(0, Snag * StepSize); normrnd(0, Snag * StepSize)];
                 PoNg(:,j) = PoNg(:,j-1) + oneStepMeasure;
                 [PoEKF(:,j), P] = MoveEKF(PoEKF(:,j-1), oneStepMeasure, P, TestRSS(i).RSSsPca(1:NapPca,j), TestRSS(i).RSSsPca(1:NapPca,j-1), eye(2), 10^rq*eye(2), AbRSSDB, RePoDB, AbRSSMinusDB, Nap, K);
                 %PoEKF(:,j) = PoEKF(:,j-1) + uEKF(:,j-1);
                 [PoFP(:,j), PFP] = FEKF(PoFP(:,j-1), oneStepMeasure, PFP, TestRSS(i).RSSsPca(1:NapPca,j), eye(2), 10^rq*eye(2), AllRSS, Nap, K);
%                  [PoEKFHigh(:,j), PHigh] = ModiEKF(PoEKFHigh(:,j-1), oneStepMeasure, oneStepMeasureBefo, PHigh,...
%                      (TestRSS(i).RSSsPca(1:NapPca,j)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), (TestRSS(i).RSSsPca(1:NapPca,j-1)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), (TestRSS(i).RSSsPca(1:NapPca,j-2)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), ...
%                      net,  eye(2), 10^rq*eye(NapPca), Abtm1ReDB, AbtDB, meanRSSNorm);
                 [PoEKFHigh(:,j), PHigh] = ModiEKF(PoEKFHigh(:,j-1), oneStepMeasure, PHigh,...
                     TestRSS(i).RSSsPca(1:NapPca,j), TestRSS(i).RSSsPca(1:NapPca,j-1),...
                     eye(2), 10^rq*eye(NapPca), Abtm1ReDB, AbtDB, Nap, K);
                 %[PoFPHigh(:,j), PFPHigh] = TEKF(PoFPHigh(:,j-1), oneStepMeasure, PFPHigh, (TestRSS(i).RSSsPca(1:NapPca,j)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), (TestRSS(i).RSSsPca(1:NapPca,j-1)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), FPnet, eye(2), 10^rq*eye(NapPca), AllRSS, meanRSSNorm);
                 [PoFPHigh(:,j), PFPHigh] = TEKF(PoFPHigh(:,j-1), oneStepMeasure, PFPHigh, TestRSS(i).RSSsPca(1:NapPca,j), eye(2), 10^rq*eye(NapPca), AllRSS, Nap, K);
%                  eu1a = [eu1a, eu1];
%                  eu2a = [eu2a, eu2];
                 %[PoFPHigh(:,j), PFPHigh] = ukf(PoFPHigh(:,j-1) + oneStepMeasure/StepSize, PFPHigh, (TestRSS(i).RSSsPca(1:NapPca,j)-RSSScaleMin)/(RSSScaleMax-RSSScaleMin), FPnet, eye(2), eye(NapPca));
                 TestPointNum = TestPointNum + 1;
                 ErrorSumNg = ErrorSumNg + norm(PoNg(:,j)-PoTrue(:,j));
                 ErrorSumEKF = ErrorSumEKF + norm(PoEKF(:,j)-PoTrue(:,j));
                 ErrorSumEKFHigh = ErrorSumEKFHigh + norm(PoEKFHigh(:,j)-PoTrue(:,j));
                 ErrorSumFP = ErrorSumFP + norm(PoFP(:,j)-PoTrue(:,j));
                 ErrorSumFPHigh = ErrorSumFPHigh + norm(PoFPHigh(:,j)-PoTrue(:,j));
                 ErrorArrayNg = [ErrorArrayNg; norm(PoNg(:,j)-PoTrue(:,j))];
                 ErrorArrayEKF = [ErrorArrayEKF; norm(PoEKF(:,j)-PoTrue(:,j))];
                 ErrorArrayEKFHigh = [ErrorArrayEKFHigh; norm(PoEKFHigh(:,j)-PoTrue(:,j))];
                 ErrorArrayFP = [ErrorArrayFP; norm(PoFP(:,j)-PoTrue(:,j))];
                 ErrorArrayFPHigh = [ErrorArrayFPHigh; norm(PoFPHigh(:,j)-PoTrue(:,j))];
                 %oneStepMeasureBefo = oneStepMeasure;
            end

%               subplot(2,2,i)
%                 
%                 hold on
%                 plot(PoTrue(1,:),PoTrue(2,:),'gs-', 'linewidth',2);
%                 hold on
%                 plot(PoNg(1,:),PoNg(2,:),'b+-','linewidth', 2)
%                 hold on
%                 plot(PoEKF(1,:),PoEKF(2,:),'mo-','linewidth', 2)
%                 hold on
%                 plot(PoFP(1,:),PoFP(2,:),'r^-','linewidth', 2)
%                 hold on
%                 plot(PoEKFHigh(1,:),PoEKFHigh(2,:),'yx-','linewidth', 2)
%                 hold on
%                 plot(PoFPHigh(1,:),PoFPHigh(2,:),'k*-','linewidth', 2)
%                 %plot(PoFPHigh(1,:),PoFPHigh(2,:),'cv-','linewidth', 2)
% %                 hold on
% %                 plot(PoTrue(1,1),PoTrue(2,1),'go-', 'linewidth',6);
% %                 hold on
% %                 plot(PoFPHigh(1,1),PoFPHigh(2,1),'ko-', 'linewidth',6);
%                 
%                 
%                 if i == 1
%                     %legend('Inert navigation','Kernel-EKF','FP');
%                     legend('Ground truth','Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');
%                     %legend('Ground truth','Inert navigation','TMS-High','Constant EKF');
%                 end
%                 pause
        end
%         EstPath(i).PoNg = PoNg;
%         EstPath(i).PoEKF = PoEKF;
    end
%     ErrorRQNg(rq+4) = ErrorSumNg / TestPointNum;
%     ErrorRQEKF(rq+4) = ErrorSumEKF / TestPointNum;
%     ErrorRQFP(rq+4) = ErrorSumFP / TestPointNum;
%     ErrorRQEKFHigh(rq+4) = ErrorSumEKFHigh / TestPointNum;
%     ErrorRQFPHigh(rq+4) = ErrorSumFPHigh / TestPointNum;
%     ErrorRQNg(K) = ErrorSumNg / TestPointNum;
%     ErrorRQEKF(K) = ErrorSumEKF / TestPointNum;
%     ErrorRQFP(K) = ErrorSumFP / TestPointNum;
%     ErrorRQEKFHigh(K) = ErrorSumEKFHigh / TestPointNum;
%     ErrorRQFPHigh(K) = ErrorSumFPHigh / TestPointNum;
%end


%save OnlineTestEmp.mat

% figure
% plot([-3:1:3],ErrorRQNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([-3:1:3],ErrorRQEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([-3:1:3],ErrorRQFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([-3:1:3],ErrorRQEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([-3:1:3],ErrorRQFPHigh / StepSize,'k*-','linewidth', 2)
% figure
% plot([1:2:21],ErrorRQNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([1:2:21],ErrorRQEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([1:2:21],ErrorRQFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([1:2:21],ErrorRQEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([1:2:21],ErrorRQFPHigh / StepSize,'k*-','linewidth', 2)

% plot([1:10],ErrorRQNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([1:10],ErrorRQEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([1:10],ErrorRQFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([1:10],ErrorRQEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([1:10],ErrorRQFPHigh / StepSize,'k*-','linewidth', 2)

% xlabel('Log Ratio'); ylabel('Locating errors(m)'); title('');
% % xlabel('K'); ylabel('Locating errors(m)'); title('');
% % %legend('Inert navigation','Kernel-EKF');
% legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

 ErrorMeanNg = ErrorSumNg / TestPointNum
 ErrorMeanEKF = ErrorSumEKF / TestPointNum
 ErrorMeanFP = ErrorSumFP / TestPointNum
 ErrorMeanEKFHigh = ErrorSumEKFHigh / TestPointNum
 ErrorMeanFPHigh = ErrorSumFPHigh / TestPointNum
 
%  figure
%  set(cdfplot(eu1a(1,:) / StepSize), 'linewidth', 1,'color', 'b', 'marker', '+');
%  hold on
%  set(cdfplot(eu1a(2,:) / StepSize), 'linewidth', 1, 'color', 'm', 'marker', 'o');
%  hold on
%  set(cdfplot(eu2a(1,:) / StepSize), 'linewidth', 1,'color', 'r', 'marker', '^');
%  hold on
%  set(cdfplot(eu2a(2,:) / StepSize), 'linewidth', 1, 'color', 'y', 'marker', 'x');
% legend('TMS-High-dx','TMS-High-dy','Constant EKF-dx','Constant EKF-dx');
% xlabel('Ratio of u'); ylabel('CDF'); title('');

figure
set(cdfplot(ErrorArrayNg / StepSize), 'linewidth', 1,'color', 'b', 'marker', '+');
hold on
set(cdfplot(ErrorArrayEKF / StepSize), 'linewidth', 1, 'color', 'm', 'marker', 'o');
hold on
set(cdfplot(ErrorArrayFP / StepSize), 'linewidth', 1, 'color', 'r', 'marker', '^');
hold on
set(cdfplot(ErrorArrayEKFHigh / StepSize), 'linewidth', 1, 'color', 'y', 'marker', 'x');
hold on
set(cdfplot(ErrorArrayFPHigh / StepSize), 'linewidth', 1, 'color', 'k', 'marker', '*');
% hold on
% set(cdfplot(ErrorArrayFPHigh / StepSize), 'linewidth', 1, 'color', 'c', 'marker', 'v');
% bar([ ErrorMeanNg/StepSize;  ErrorMeanEKFHigh/StepSize; ErrorMeanFPHigh/StepSize])
% set(PFCdf, 'linewidth', 2, 'color', 'k', 'marker', 'x');
% set(sdtwSLCdf, 'linewidth', 2, 'color', 'y', 'marker', '^');
% set(CosineCdf, 'linewidth', 2, 'color', 'm', 'marker', '*');
% legend('Inert navigation','Kernel-EKF','FP','Kernel-EKF-High','FP-High');
% legend('Inert navigation','TMS-High','Constant EKF');
legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');
xlabel('Locating errors(m)'); ylabel('CDF'); title('');

% set(gca,'FontSize',14)
% load OfflineTrain.mat
% [mx, my] = meshgrid(linspace(0,100,100),linspace(0,100,100));
% for i = 1 : size(mx,1)
%     for j = 1 : size(mx,2)
%         d = pdist2(ApLoc(6,:), [mx(i,j), my(i,j)]);
%         mz(i,j) = St-S0-10*2*log(d/1) - 200;
%     end
% end
% surf(mx,my,mz);
