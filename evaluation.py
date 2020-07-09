#%%
import numpy as np
from sklearn import metrics

num_srv = 32
#%%
# load test classes:
Y_test = np.load('experiments//y_test.npy')

# number of experiments:
num_existing_links = len(Y_test)

# load result classes:
Y_pred = np.zeros(num_existing_links,bool)
Y_exp_pred = np.zeros(num_existing_links,bool)

for i in range(num_srv):
    Y_exp_pred = np.load('experiments//y_pred_'+str(i+1)+'.npy')
    Y_pred+=Y_exp_pred
#%%
# cost-based validation:
validation_cost = np.zeros(len(Y_test),float)
validation_cost+=1
if np.abs(len(Y_test)-np.sum(Y_test)) != 0:
    neg_weight = len(Y_test)/float(2)/float(np.abs(len(Y_test)-np.sum(Y_test)))
    for i in range(len(Y_test)):
        if Y_test[i] == False:
            validation_cost[i] = neg_weight

# quality measures:
false_positive,true_positive,thresholds = metrics.roc_curve(Y_test,Y_pred)
accuracy = metrics.accuracy_score(Y_test, Y_pred)
acc_weighted = metrics.accuracy_score(Y_test, Y_pred, sample_weight=validation_cost)
auc_score = metrics.auc(false_positive,true_positive)
precision = metrics.precision_score(Y_test,Y_pred, average='binary')
prec_weighted = metrics.precision_score(Y_test,Y_pred, average='binary', sample_weight=validation_cost)
recall = metrics.recall_score(Y_test,Y_pred,average='binary')
rec_weighted = metrics.recall_score(Y_test,Y_pred,average='binary', sample_weight=validation_cost)

print 'Accuracy: ',accuracy
print 'Accuracy (weighted): ',acc_weighted
print 'AUC score: ',auc_score
print 'Precision: ',precision
print 'Precision (weighted): ',prec_weighted
print 'Recall: ',recall
print 'Recall (weighted): ',rec_weighted
#%%
counts = []
times = []

for i in range(num_srv):
    f = open('experiments//y_pred_'+str(i+1)+'.txt','r')
    for line in f:
        try:
            if 'PageRank count' in line:
                data = line.split(':')
                counts.append(np.int(data[1]))
            elif 'Time spent, sec' in line:
                data = line.split(':')
                times.append(np.float(data[1]))
        except ValueError:
            print "Invalid input:", line
    f.close()
    
np.min(counts)
np.max(counts)
np.mean(counts)
np.median(counts)

np.min(times)
np.max(times)
np.mean(times)
np.median(times)
#%%
