#!/usr/bin/env python

"""Machine Learning Classifiers

Functions for running and evaluating ML classifiers

"""

from time import time
from collections import defaultdict
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.utils.extmath import density
from sklearn.feature_selection import SelectFromModel
from sklearn import metrics
from sklearn.linear_model import RidgeClassifier
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import Perceptron
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import RandomForestClassifier

import matplotlib.pyplot as plt

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 160 else s[:157] + "..."


def classifiers():
    """List of sklearn classifiers to test


    """

    clfs = [
        (RidgeClassifier(tol=1e-2, solver="sag"),'RidgeClassifier'),
        (Perceptron(n_iter=50),'Perceptron'),
        (PassiveAggressiveClassifier(n_iter=50),'PassiveAggressiveClassifier'),
        (KNeighborsClassifier(n_neighbors=10),'KNeighborsClassifier'),
        (RandomForestClassifier(n_estimators=100),'RandomForestClassifier'),
        (LinearSVC(penalty="l1", dual=False, tol=1e-3),'LinearSVC_l1'),
        (LinearSVC(penalty="l2", dual=False, tol=1e-3),'LinearSVC_l2'),
        (SGDClassifier(alpha=.0001, n_iter=50, penalty="l1"),'SGDClassifier_l1'),
        (SGDClassifier(alpha=.0001, n_iter=50, penalty="l2"),'SGDClassifier_l2'),
        (SGDClassifier(alpha=.0001, n_iter=50, penalty="elasticnet"),'SGDClassifier_elastinet'),
        (NearestCentroid(),'NearestCentroid'),
        (MultinomialNB(alpha=.01),'MultinomialNB'),
        (BernoulliNB(alpha=.01),'BernoulliNB'),
        (Pipeline([
            ('feature_selection', SelectFromModel(LinearSVC(penalty="l1", dual=False, tol=1e-3))),
            ('classification', LinearSVC(penalty="l2"))
        ]),'LinearSVC_modelselection')
    ]

    return clfs


def benchmark(classifier, X_train, y_train, X_test, y_test, logger, feature_names, target_names, opts):
    """Benchmark sklearn classifiers


    """

    logger.info('_' * 80)
    logger.info("Training: ")
    logger.info(classifier)
    t0 = time()
    classifier.fit(X_train, y_train)
    train_time = time() - t0
    logger.info("train time: %0.3fs" % train_time)

    t0 = time()
    pred = classifier.predict(X_test)
    test_time = time() - t0
    logger.info("test time:  %0.3fs" % test_time)

    score = metrics.accuracy_score(y_test, pred)
    logger.info("accuracy:   %0.3f" % score)

    if hasattr(classifier, 'coef_'):
        logger.info("dimensionality: %d" % classifier.coef_.shape[1])
        logger.info("density: %f" % density(classifier.coef_))

        if opts.info_top10 and feature_names is not None:
            logger.info("top 10 regions per class:")
            for i, label in enumerate(target_names[1:]):
                top10 = np.argsort(classifier.coef_[i])[-10:]
                logger.info(trim("%s: %s" % (label, " ".join(feature_names[top10]))))
        logger.info("\n")

    if opts.info_report:
        logger.info("classification report:")
        logger.info(metrics.classification_report(y_test, pred,
                                            target_names=target_names))
    if opts.info_cm:
        logger.info("confusion matrix:")
        logger.info(metrics.confusion_matrix(y_test, pred))

    logger.info("\n")
    logger.info('=' * 80)
    classifier_descr = str(classifier).split('(')[0]
    return classifier_descr, score, train_time, test_time


def barplot(results):
    """Plot runtime and success rate for each classifier

    results: list of tuples containing classifier name, success rate, training time, prediction time

    """

    indices = np.arange(len(results))
    results = [[x[i] for x in results] for i in range(4)]

    clf_names, score, training_time, test_time = results
    training_time = np.array(training_time) / np.max(training_time)
    test_time = np.array(test_time) / np.max(test_time)

    plt.figure(figsize=(12, 8))
    plt.title("Score")
    plt.barh(indices, score, .2, label="score", color='navy')
    plt.barh(indices + .3, training_time, .2, label="training time",
             color='c')
    plt.barh(indices + .6, test_time, .2, label="test time", color='darkorange')
    plt.yticks(())
    plt.legend(loc='best')
    plt.subplots_adjust(left=.25)
    plt.subplots_adjust(top=.95)
    plt.subplots_adjust(bottom=.05)

    for i, c in zip(indices, clf_names):
        plt.text(-.3, i, c)

    plt.show()


def run(pg, amr, logger, feature_names, target_names, drug_names, opts):
    """Run 10 tests for each amr phenotype

    """

    n = 10
    ndrug = amr.shape[1]
    clfs = classifiers()
    results = {}

    for d in range(ndrug):
        validrows = ~np.isnan(amr[:,d])

        X = pg[validrows,:]
        y = amr[validrows,d]

        drug = drug_names[d]

        results[drug] = {}
        for c,cname in clfs:
            results[drug][cname] = {'score': [], 'train_time': [], 'test_time': [] }

        for i in range(n):
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i+3)

            for c,cname in clfs:
                _, score, train_time, test_time = benchmark(c, X_train, y_train, X_test, y_test, logger, feature_names, target_names, opts)
                results[drug][cname]['score'].append(score)
                results[drug][cname]['train_time'].append(train_time)
                results[drug][cname]['test_time'].append(test_time)


    return results



def summaryplot(results, logger=None):
    """Plot mean success rate for each classifier per drug

    

    """
    x = results.keys()

    y = defaultdict(list)
    errs = defaultdict(list)
    for d in x:
        for c in results[d]:
            arr = np.array(results[d][c]['score'])

            y[c].append(arr.mean())
            errs[c].append(arr.std())

    if logger:
        logger.info("Accuracy Summary Report: ")
        

    n = len(x)
    ind = np.arange(n)
    
    clfrs = list(y.keys())
    blocksize = 7
    width = 1/(blocksize+1)
    subplots = [clfrs[i:i + blocksize] for i in range(0, len(clfrs), blocksize)]

    for l in subplots:

        fig, ax = plt.subplots()
        i = 0
        bars = []
        labels = []
        for c in l:
            rects = ax.bar(ind+width*i, y[c], width, yerr=errs[c])
            bars.append(rects[0])
            labels.append(c)
            i=i+1

            if logger:
                info_list = ['{}: {}'.format(s[0],np.round(s[1],4)) for s in zip(x,y[c])]
                logger.info('_' * 80)
                logger.info('{} Values\n{}'.format(c, ', '.join(info_list)))
                logger.info('=' * 80)

        ax.legend(bars, labels, loc=('lower left'))
        ax.set_ylabel('Accuracy Scores')
        ax.set_title('Accuracy for different ML Classifiers by Antimicrobial Drug')
        ax.set_yticks(np.arange(0,1.05,0.05))
        #ax.set_xticks(ind + width*n / 2)
        ax.set_xticks(ind)
        ax.set_xticklabels(x)
        ax.grid(linestyle='--')
        xlabels = ax.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=10)

    plt.show()











