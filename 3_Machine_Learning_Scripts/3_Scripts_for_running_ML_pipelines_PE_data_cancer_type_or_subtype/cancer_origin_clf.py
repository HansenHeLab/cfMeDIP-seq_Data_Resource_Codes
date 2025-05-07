"""
File:      cancer_origin_clf.py
Author:    Yuanchang Fang
Date:      March 2025
Purpose:
    1. Load engineered fragmentomic and methylation feature matrices for pan-cancer classification.
    2. Harmonize sample IDs against metadata (filtering to primary PE cohort and removing normals/HCC).
    3. Build and evaluate multiple classifiers (LR, SVM, KNN, RF, XGBoost, GBM, LDA) within a nested
       cross‑validation framework (Stratified K‑Fold outer loop, inner GridSearch for hyperparameter tuning).
    4. Automatically adapt CV folds based on class imbalance, and compute ROC AUC plus sensitivities at
       specified specificity thresholds (80%, 85%, 90%, 95%, 99%).
    5. Export per‑cancer & per‑feature summary tables (.sum.csv) and detailed ROC curve data (.roc.csv) for all runs.

Description:
    This script performs **one‑versus‑all binary classification** for each individual cancer type
    (e.g., “Breast Cancer”, “Lung Cancer”, “AML”, etc.) using multiple fragmentomic and 
    methylation‑derived feature sets. For each target cancer, all other samples are treated as 
    the negative class, and classifiers are trained/evaluated in a nested cross‑validation 
    framework that adapts fold counts based on the prevalence of the target class.

Functionality:
    1. **Data Loading & Assembly**  
       - Reads one or more feature matrices (motif, methylation, DELFI ratios, insert‑size, 
         nucleosome peaks) from CSV files, inner‑merging them on sample IDs.  
       - Normalizes sample ID types and aligns with metadata.

    2. **Metadata Filtering & Labeling**  
       - Loads master metadata (cancer type, project, sequencing type).  
       - Filters out healthy (“Normal”) and unwanted cohorts (e.g., TCGE‑CFMe‑HCC), retaining 
         primary PE samples.  
       - Creates binary labels `y`:  
           • `1` = samples of the specified cancer type  
           • `0` = all other samples (one‑vs‑all)

    3. **Nested Cross‑Validation & Model Selection**  
       - **Outer loop**: Stratified K‑Fold (fold count chosen by class balance) to hold out test sets.  
       - **Inner loop**: Stratified K‑Fold GridSearch to jointly tune PCA component count and 
         classifier hyperparameters.  
       - **Classifiers** supported:
           • Logistic Regression (L1)  
           • Support Vector Machine (poly, rbf, linear)  
           • K‑Nearest Neighbors  
           • Random Forest  
           • XGBoost  
           • Gradient Boosting Machine  
           • Linear Discriminant Analysis  
       - Automatically handles class imbalance via class/sample weighting.

    4. **Performance Metrics**  
       - **Primary**: ROC AUC on held‑out test folds (one‑vs‑all).  
       - **Secondary**: Sensitivity at fixed specificity thresholds (80%, 85%, 90%, 95%, 99%).  
       - Aggregates per‑fold predictions to build detailed ROC data.

    5. **Outputs**  
       - **Summary Tables** (`.sum.csv`): one row per fold/run, including:
           • `feature`, `cancer`, `model`, `iteration`, `run`  
           • Best inner‑loop PCA components & hyperparameters  
           • Validation AUC, test AUC  
           • Sensitivities at each specificity threshold  
           • `target_num` (positive) & `nontarget_num` (negative) counts  
       - **ROC Curve Data** (`.roc.csv`): concatenated `fpr` & `tpr` values across test folds,
         enabling per‑cancer one‑vs‑all ROC plotting downstream.

Usage:
    $ python run_classification.py <feature> "<cancer_type>"
    e.g.:
      python run_classification.py motif "Breast Cancer"

Environment:
    - Python 3.x  
    - Packages: numpy, pandas, scikit‑learn, xgboost, matplotlib
"""


import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, LeaveOneOut
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, confusion_matrix, roc_curve, auc
from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight
from sklearn.decomposition import PCA
import xgboost as xgb
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sys

import os
print("Total CPUs:", os.cpu_count())

seed = 123

def get_data(data_sets):
    if len(data_sets) == 1:
        X = pd.read_csv(data_sets[0], index_col = 0)
    else:
        first = True
        for data in data_sets:
            data = pd.read_csv(data, index_col = 0)
            if first:
                X = data.copy()
                first = False
            else:
                # inner concat to ensure no NAs
                X = pd.concat([X, data], join = 'inner', axis = 1)
    # index is a mix of different dtypes, need to convert to string
    X.index = X.index.astype(str)
    return X

# def train_with_loo(X, y, clf_name):

#     out_dict = {'best_params': [], 
#                 'best_val_accuracy': [], 
#                 'iteration': [],
#                 'run': []} 
    
#     run_arr = np.array([])
#     fpr_arr = np.array([])
#     tpr_arr = np.array([])

#     pipe = Pipeline([
#         ('scaler', StandardScaler()),
#         ('pca', PCA()),
#         (clf_name, classifiers[clf_name])
#         ])

#     run = '.'.join([str(x) for x in [cancer, feature, clf_name, 1]])

#     loo = LeaveOneOut()
#     grid_search = GridSearchCV(estimator = pipe, param_grid = {**grids['pca'], **grids[clf_name]}, cv = loo, n_jobs = -1, verbose = 1, return_train_score = False, scoring = 'accuracy')

#     if clf_name in ['xgb','gbm']:
#         # calculate sample weight
#         sample_w = compute_sample_weight(class_weight = 'balanced', y = y)
#         grid_search.fit(X, y, **{f'{clf_name}__sample_weight': sample_w})
#     else:
#         grid_search.fit(X, y)

#     best_params = grid_search.best_params_
#     best_score = grid_search.best_score_
#     best_clf = grid_search.best_estimator_

#     out_dict['best_params'].append(best_params)
#     out_dict['best_val_accuracy'].append(best_score)
#     out_dict['iteration'].append(1)
#     out_dict['run'].append(run)

#     # get roc on the whole date
#     y_pred_prob = best_clf.predict_proba(X)[:, 1]
#     fpr, tpr, _ = roc_curve(y, y_pred_prob)
#     run_arr = np.append(run_arr, [run] * len(fpr))
#     fpr_arr = np.append(fpr_arr, fpr)
#     tpr_arr = np.append(tpr_arr, tpr)
#     roc_data = {'run': run_arr, 'fpr': fpr_arr, 'tpr': tpr_arr} 

#     return out_dict, roc_data

def get_sensitivity_at_specificity(fpr, tpr, desired_specificity):
    specificity = 1 - fpr 
    # nonzero returns the index of the non-zero elements and it's a tuple
    # specificity is decreasing, so we want the last element that is greater than or equal to the desired specificity
    idx = (specificity >= desired_specificity).nonzero()[0][-1]
    target_sensitivity = tpr[idx]
    return target_sensitivity

def train_with_nkf(X, y, clf_name, n, k):

    out_dict = {'best_params': [], 
                'best_nPCs': [],
                'best_val_auc': [], 
                'test_auc': [], 
                'test_sensitivity_80specificity': [],
                'test_sensitivity_85specificity': [],
                'test_sensitivity_90specificity': [],
                'test_sensitivity_95specificity': [],
                'test_sensitivity_99specificity': [],
                'iteration': [],
                'run': []} 
    
    run_arr = np.array([])
    fpr_arr = np.array([])
    tpr_arr = np.array([])

    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA()),
        (clf_name, classifiers[clf_name])
        ])

    cv_outer = StratifiedKFold(n_splits = n, shuffle = True, random_state = seed)
    iteration = 0
    for train_outer_index, test_index in cv_outer.split(X, y):
        iteration += 1
        run = '.'.join([str(x) for x in [cancer, feature, clf_name, iteration]])

        X_outer_train, y_outer_train = X[train_outer_index, :], y[train_outer_index]
        X_test, y_test = X[test_index, :], y[test_index]

        cv_inner = StratifiedKFold(n_splits = k, shuffle = True, random_state = seed)
        grid_search = GridSearchCV(estimator = pipe, param_grid = {**grids['pca'], **grids[clf_name]}, cv = cv_inner, n_jobs = -1, verbose = 1, return_train_score = False, scoring = 'roc_auc')

        if clf_name in ['xgb','gbm']:
            # calculate sample weight
            sample_w = compute_sample_weight(class_weight = 'balanced', y = y_outer_train)
            grid_search.fit(X_outer_train, y_outer_train, **{f'{clf_name}__sample_weight': sample_w})
        else:
            grid_search.fit(X_outer_train, y_outer_train)

        best_inner_params = grid_search.best_params_
        best_inner_score = grid_search.best_score_
        best_clf = grid_search.best_estimator_

        # get number of PCs
        nPCs = best_clf.named_steps['pca'].n_components_  

        # now test for this iteration
        y_pred_prob = best_clf.predict_proba(X_test)[:, 1]
        test_score = roc_auc_score(y_test, y_pred_prob)

        out_dict['best_params'].append(best_inner_params)
        out_dict['best_nPCs'].append(nPCs)
        out_dict['best_val_auc'].append(best_inner_score)
        out_dict['test_auc'].append(test_score)
        out_dict['iteration'].append(iteration)
        out_dict['run'].append(run)

        fpr, tpr, _ = roc_curve(y_test, y_pred_prob)

        sen_80spec = get_sensitivity_at_specificity(fpr, tpr, 0.8)
        sen_85spec = get_sensitivity_at_specificity(fpr, tpr, 0.85)
        sen_90spec = get_sensitivity_at_specificity(fpr, tpr, 0.9)
        sen_95spec = get_sensitivity_at_specificity(fpr, tpr, 0.95)
        sen_99spec = get_sensitivity_at_specificity(fpr, tpr, 0.99)

        out_dict['test_sensitivity_80specificity'].append(sen_80spec)
        out_dict['test_sensitivity_85specificity'].append(sen_85spec)
        out_dict['test_sensitivity_90specificity'].append(sen_90spec)
        out_dict['test_sensitivity_95specificity'].append(sen_95spec)
        out_dict['test_sensitivity_99specificity'].append(sen_99spec)

        run_arr = np.append(run_arr, [run] * len(fpr))
        fpr_arr = np.append(fpr_arr, fpr)
        tpr_arr = np.append(tpr_arr, tpr)
        roc_data = {'run': run_arr, 'fpr': fpr_arr, 'tpr': tpr_arr} 

    return out_dict, roc_data

classifiers = {
    'lr': LogisticRegression(class_weight = 'balanced', random_state = seed),
    'svm': SVC(class_weight = 'balanced', probability = True, random_state = seed),
    'knn': KNeighborsClassifier(), # doesn't support class_weight and random state
    'rf': RandomForestClassifier(class_weight = 'balanced', random_state = seed),
    'xgb': xgb.XGBClassifier(random_state = seed), # add class weight in fit
    'gbm': GradientBoostingClassifier(random_state = seed), # add class weight in fit
    'lda': LinearDiscriminantAnalysis() # doesn't support class_weight and random state
}

pca_grid = {
    'pca__n_components': [0.6, 0.7, 0.8, 0.9, 0.95]
}

lr_grid = {
    'lr__penalty': ['l1'],
    'lr__solver': ['saga'],
    'lr__C': np.linspace(0.0001, 1, 10)
    }

svm_grid = {
    'svm__kernel': ['poly', 'rbf', 'linear'],
    'svm__C': [0.01, 0.1, 1, 10, 100]
    }

knn_grid = {
    'knn__n_neighbors': [5, 10, 50, 100],
    'knn__weights': ['uniform', 'distance']
    }

rf_grid = {
    'rf__max_depth': [5, 20, None],
    'rf__n_estimators': [50, 100, 500]
    }


xgb_grid = {
        'xgb__gamma': [0, 1],
        'xgb__max_depth': [3, 6],
        'xgb__n_estimators': [100, 200],
        'xgb__learning_rate': [0.01, 0.1],
        'xgb__colsample_bytree': [0.5, 0.8],
        'xgb__min_child_weight': [1, 3],
        'xgb__subsample': [0.7, 1]
        }

gbm_grid = {
    'gbm__max_depth': [2, 3, 4],          
    'gbm__n_estimators': [100, 150, 200], 
    'gbm__learning_rate': [0.01, 0.1],    
    'gbm__min_samples_leaf': [3, 5]
}

lda_grid = {
    'lda__solver': ['svd', 'lsqr', 'eigen']
}

grids = {
    'pca': pca_grid,
    'lr': lr_grid,
    'svm': svm_grid,
    'knn': knn_grid,
    'rf': rf_grid,
    'xgb': xgb_grid,
    'gbm': gbm_grid,
    'lda': lda_grid
    }

main_dir = '/.mounts/labs/PCSI/users/yfang/misc_projects/cfMEDIP/ML'

meta_data = pd.read_csv(f'{main_dir}/data/metadata_df.csv')
# remove normal samples and the HCC project, and use PE samples only
meta_data = meta_data[(meta_data['cancer_type'] != 'Normal') & (meta_data['project_id'] != 'TCGE-CFMe-HCC') & (meta_data['type'] == 'PE')]
meta_data['sample_id'] = meta_data['sample_id'].str.replace('_dedup', '')
meta_data.set_index('sample_id', inplace = True)
# index is a mix of different dtypes, need to convert to string
meta_data.index = meta_data.index.astype(str)
meta_data_samples = set(meta_data.index)

# input data
motif_data = f'{main_dir}/data/motif_merged.csv'
methyl_data = f'{main_dir}/data/methyl_merged.csv'
delfi_data = f'{main_dir}/data/delfi_merged.csv'
insert_size_data = f'{main_dir}/data/insert_size_merged.csv'
ns_peaks_data = f'{main_dir}/data/ns_peaks_merged.csv'

# nkf has best_val_auc, test_auc, test_sensitivity, test_specificity
sum_tab = pd.DataFrame(columns = ['feature', 'cancer', 'model', 'iteration', 'run', 'target_num', 'nontarget_num', 'best_params', 'best_nPCs', 'best_val_auc',
                                  'test_sensitivity_80specificity', 'test_sensitivity_85specificity', 'test_sensitivity_90specificity',
                                  'test_sensitivity_95specificity', 'test_sensitivity_99specificity', 'test_auc'])
roc_data_sum = pd.DataFrame(columns = ['run', 'fpr', 'tpr'])

"""
Features:
'motif', 'methyl', 'delfi', 'insert_size', 'ns_peaks', 'motif+methyl', 'motif+methyl+delfi', 'delfi+methyl', 'all'
"""
feature = sys.argv[1]

if feature in ['motif', 'methyl', 'delfi', 'insert_size', 'ns_peaks']:
    X = get_data([eval(f'{feature}_data')])
elif feature == 'motif+methyl':
    X = get_data([motif_data, methyl_data])
elif feature == 'motif+methyl+delfi':
    X = get_data([motif_data, methyl_data, delfi_data])
elif feature == 'delfi+methyl':
    X = get_data([delfi_data, methyl_data])
elif feature == 'all':
    X = get_data([motif_data, methyl_data, delfi_data, insert_size_data, ns_peaks_data])

# find samples that are in both meta data and the feature table
feature_table_samples = set(X.index)
final_samples = list(meta_data_samples & feature_table_samples)

# ensure the order is the same
X = X.loc[final_samples]
X = np.array(X)

meta_data = meta_data.loc[final_samples]

"""
All cancer types:
'Brain Cancer', 'Head and Neck Cancer', 'AML', 'Prostate Cancer', 'Lung Cancer', 'Eye Cancer', 'LFS Survivor', 'LFS Previvor',
'LFS Positive', 'Breast Cancer', 'Ovarian Cancer', 'Melanoma', 'Mixed Cancer'
"""
cancer = sys.argv[2]

y = np.array(meta_data['cancer_type'])
y = np.where(y == cancer, 1, 0)

target_num = y.sum()
nontarget_num = len(y) - target_num

for clf in classifiers.keys():
    
    print('#'*50)
    print(f'feature = {feature}, cancer = {cancer}, clf = {clf}')

    smaller_n = min(y.sum(), len(y) - y.sum())
    if smaller_n < 10:
        # print('LOOCV')
        # out_res = train_with_loo(X, y, clf)
        print('NKF 5*3')
        out_res = train_with_nkf(X, y, clf, n = 5, k = 3)
    elif 10 <= smaller_n <= 100:
        print('NKF 5*5')
        out_res = train_with_nkf(X, y, clf, n = 5, k = 5)
    elif 100 < smaller_n <= 250:
        print('NKF 10*5')
        out_res = train_with_nkf(X, y, clf, n = 10, k = 5)
    else:
        print('NKF 10*10')
        out_res = train_with_nkf(X, y, clf, n = 10, k = 10)
        

    sub_out = pd.DataFrame(out_res[0])
    sub_out[['feature', 'cancer', 'model', 'target_num', 'nontarget_num']] = [feature, cancer, clf, target_num, nontarget_num]
    sum_tab = pd.concat([sum_tab, sub_out], axis = 0)

    roc_data = out_res[1]
    roc_data_sum = pd.concat([roc_data_sum, pd.DataFrame(roc_data)], axis = 0)

    print('$'*100)
    print(len(roc_data_sum))

sum_tab.to_csv(f'{main_dir}/results_cancers/{cancer}+{feature}.sum.csv', index = False)
roc_data_sum.to_csv(f'{main_dir}/results_cancers/{cancer}+{feature}.roc.csv', index = False)