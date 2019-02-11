global Snag;
global ErrorMeanNg; 
global ErrorMeanEKF;
global ErrorMeanFP;
global ErrorMeanEKFHigh;
global ErrorMeanFPHigh;
global ErrorMeanArrayNg;
global ErrorMeanArrayEKF;
global ErrorMeanArrayFP;
global ErrorMeanArrayEKFHigh;
global ErrorMeanArrayFPHigh;
global ErrorNgTt;
global ErrorEKFTt;
global ErrorFPTt;
global Svar;
global PathNum;
global Nap;
global Restart;
global SnagTrain;
global SnagTest;
global ErrorArrayNN;
global ErrorNN;
global PathTrain;
global PathTest;
global ApLoc;
global NNnum;
global NapPca;
global tt;
global mapping;
global emp;
global ReScale;
global trainrandp;
global testrandp;

ErrorMeanArrayNg = [];
ErrorMeanArrayEKF = [];
ErrorMeanArrayFP = [];
ErrorMeanArrayEKFHigh = [];
ErrorMeanArrayFPHigh = [];
ErrorArrayNN = [];
%Restart = 1;

% for Snag = 0: 0.1: 1
%     Snag
%     if Snag == 0
%         Restart = 1;
%      else
%         Restart = 1;
%      end
%     OfflineTrain;
%     %OfflineNN;
%     %OfflineNNFP;
%     OfflineNNMove;
%     OnlineTest;
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
%     ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
%     ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
% end
% figure
% plot([0:0.1:1],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([0:0.1:1],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([0:0.1:1],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([0:0.1:1],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([0:0.1:1],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
% xlabel('V_u'); ylabel('Locating errors(m)'); title('');
% legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

% for ReScale = 0: 1: 10%: 30
%     if ReScale == 0
%         Restart = 1;
%     else
%         Restart = 0;
%     end
%     OfflineTrain;
%     %OfflineNN;
%     %OfflineNNFP;
%     OfflineNNMove;
%     OnlineTest;
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
%     ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
%     ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
%     ErrorArrayNN = [ErrorArrayNN, ErrorNN]
% end
% 
% figure
% plot([0:1:10],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([0:1:10],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([0:1:10],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([0:1:10],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([0:1:10],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
% xlabel('V_z'); ylabel('Locating errors(m)'); title('');
% %legend('Inert navigation','Kernel-EKF','Fingerprint');
% legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');


% for Svar = 0: 10: 100%: 30
%     if Svar == 0
%         Restart = 1;
%     else
%         Restart = 0;
%     end
%     OfflineTrain;
%     %OfflineNN;
%     %OfflineNNFP;
%     OfflineNNMove;
%     OnlineTest;
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
%     ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
%     ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
%     ErrorArrayNN = [ErrorArrayNN, ErrorNN]
% end
% 
% figure
% plot([0:10:100],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([0:10:100],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([0:10:100],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([0:10:100],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([0:10:100],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
% xlabel('V_z'); ylabel('Locating errors(m)'); title('');
% %legend('Inert navigation','Kernel-EKF','Fingerprint');
% legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

for emp = 0: 0.1: 1%: 30
    if emp == 0
        Restart = 1;
    else
        Restart = 0;
    end
    OfflineTrain;
%     OfflineNN;
%     OfflineNNFP;
    OfflineNNMove;
    OnlineTest;
    ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
    ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
    ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
    ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
    ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
    ErrorArrayNN = [ErrorArrayNN, ErrorNN]
end

figure
plot([0:0.1:1],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
hold on
plot([0:0.1:1],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
hold on
plot([0:0.1:1],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
hold on
plot([0:0.1:1],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
hold on
plot([0:0.1:1],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
xlabel('Empty Ratio'); ylabel('Locating errors(m)'); title('');
legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

% ErrorNgTt = zeros(1,7);
% ErrorEKFTt = zeros(1,7);
% ErrorFPTt = zeros(1,7);
% tt = 0;
% for tt = 1 : 1%100
% ErrorMeanArrayNg = [];
% ErrorMeanArrayEKF = [];
% ErrorMeanArrayFP = [];
% ErrorArrayNN = [];
% for PathNum = 101: -10: 1
%     PathNum
%     if PathNum == 101
%           Restart = 1;
%     else
%           Restart = 0;
%      end
%     OfflineTrain;
%     %OfflineNN;
%     %OfflineNNFP;
%     OfflineNNMove;
%     OnlineTest;
%     PathNum
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
%     ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
%     ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
%     %ErrorArrayNN = [ErrorArrayNN, ErrorNN]
% end
% %     ErrorNgTt  = ErrorMeanArrayNg + ErrorNgTt
% %     ErrorEKFTt  = ErrorMeanArrayEKF + ErrorEKFTt
% %     ErrorFPTt  = ErrorMeanArrayFP + ErrorFPTt
% % end
% 
% figure
% plot([101:-10:1],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([101:-10:1],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([101:-10:1],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([101:-10:1],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([101:-10:1],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
% xlabel('Training Data Size'); ylabel('Locating errors(m)'); title('');
% %legend('Inert navigation','Kernel-EKF','Fingerprint');
% legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

% for NapPca = 1 : 5 : 26
%     if NapPca == 1
%        Restart = 1;
%     else
%         Restart = 0;
%     end
%     OfflineTrain;
%     OfflineNN;
%     OfflineNNFP;
%     OnlineTest;
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
% end
% figure
% plot([1:5:26],ErrorMeanArrayNg / StepSize,'bo-','linewidth', 2)
% hold on
% plot([1:5:26],ErrorMeanArrayEKF / StepSize,'mx-','linewidth', 2)
% hold on
% plot([1:5:26],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% xlabel('The number of AP'); ylabel('Locating errors(m)'); title('');
% legend('Inert navigation','Kernel-EKF','Fingerprint');

% for Nap = 1 : 5 : 51
%     Nap
%     if Nap == 1
%        Restart = 1;
%     else
%         Restart = 0;
%     end
%     OfflineTrain;
% %     OfflineNN;
% %     OfflineNNFP;
%     OfflineNNMove;
%     OnlineTest;
%     ErrorMeanArrayNg = [ErrorMeanArrayNg, ErrorMeanNg]
%     ErrorMeanArrayEKF = [ErrorMeanArrayEKF, ErrorMeanEKF]
%     ErrorMeanArrayFP = [ErrorMeanArrayFP, ErrorMeanFP]
%     ErrorMeanArrayEKFHigh = [ErrorMeanArrayEKFHigh, ErrorMeanEKFHigh]
%     ErrorMeanArrayFPHigh = [ErrorMeanArrayFPHigh, ErrorMeanFPHigh]
% end
% 
% plot([1:5:51],ErrorMeanArrayNg / StepSize,'b+-','linewidth', 2)
% hold on
% plot([1:5:51],ErrorMeanArrayEKF / StepSize,'mo-','linewidth', 2)
% hold on
% plot([1:5:51],ErrorMeanArrayFP / StepSize,'r^-','linewidth', 2)
% hold on
% plot([1:5:51],ErrorMeanArrayEKFHigh / StepSize,'yx-','linewidth', 2)
% hold on
% plot([1:5:51],ErrorMeanArrayFPHigh / StepSize,'k*-','linewidth', 2)
% xlabel('The number of AP'); ylabel('Locating errors(m)'); title('');
% %legend('Inert navigation','Kernel-EKF','Fingerprint');
%legend('Inert navigation','TMS-Low','FP-Low','TMS-High','FP-High');

% save OnlineTestNap.mat