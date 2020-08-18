tpr=[1.0000000 0.9622919 0.0000000]
fpr=[1.0000000 0.1189262 0.0000000]

plot(fpr,tpr,'-ow')

set(groot,'defaultLineLineWidth',1.0)

hold on

plot([0 1],[0 1],'--b')

set(gca,'FontSize',10,'FontWeight','bold')

xlabel('False Positive Rate','FontSize',12,'FontWeight','bold')
ylabel('True Positive Rate','FontSize',12,'FontWeight','bold')