test=[93.29 94.34 93.75;	
	94.02 95.95 96.46;	
	 92.64 94.62	95.36;
		91.59	93.89	95.13;
		91.13	93.89	95.13]

train=[	92.2	93.7	93;
		94.8	96.7	97.7;
		97	99.6	99.9;
		99.4	100	100;
		100	100	100]
		
C = [0.001,0.01,0.1,1,10];

figure();
semilogx(C,test(:,1)/100,C,train(:,1)/100);
legend('Test acc.','Train acc.');
title('Prediction accuracy m=1');
ylabel('Acc.');
xlabel('C');
ylim([0.91,1]);

figure();
semilogx(C,test(:,2)/100,C,train(:,2)/100);
legend('Test acc.','Train acc.');
title('Prediction accuracy m=2');
ylabel('Acc.');
xlabel('C');
ylim([0.91,1]);

figure();
semilogx(C,test(:,3)/100,C,train(:,3)/100);
legend('Test acc.','Train acc.');
title('Prediction accuracy m=3');
ylabel('Acc.');
xlabel('C');
ylim([0.91,1]);

figure();
semilogx(C,test(:,1)/100,C,test(:,2)/100,C,test(:,3)/100);
legend('m=1','m=2','m=3');
title('Prediction accuracy on testing data set')
ylabel('Acc.');
xlabel('C');
