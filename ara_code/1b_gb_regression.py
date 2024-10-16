"""
Gradient Boosting Regression

# Required Inputs
    -X      Path to feature matrix file
    -y_name Column name of label in X matrix
    -test   File containing list of test instances
    -save   Path to save output files
    -prefix Prefix of output file names
    
    # Optional
    -Y      Path to label matrix file, if label not in X matrix
    -tag    Feature types/identifier for output file naming
    -fold   k folds for Cross-Validation (default is 5)
    -n      Number of CV repetitions (default is 10)
    -feat   File containing features (from X) to include in model
    -plot   Plot feature importances and predictions (default is t)

# Outputs for each training repetition (prefixed with <prefix>_)
    _lm_test_rep_*.pdf        Regression plot of predicted and actual test labels
    _model_rep_*.save         Gradient boosting model
    _top20_rep_*.pdf          Plot of top 20 features' importance scores

# Summary outputs (prefixed with <prefix>_)
    _imp.csv                  Feature importance scores
    _cv_results.csv           Cross-validation results (various metrics)
    _test_results.csv         Evaluation results (various metrics)
    RESULTS_gradboost.txt       Aggregated results (various metrics)
"""
__author__ = "Kenia Segura Abá"

from configparser import ExtendedInterpolation
import sys, os, argparse
import time
import random
import pickle
import datatable as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from hyperopt import hp, fmin, tpe, Trials
from hyperopt.pyll.base import scope
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score, explained_variance_score
from sklearn.model_selection import KFold, cross_validate
from sklearn.model_selection import cross_val_predict


def hyperopt_objective(params):
    """
    Create the hyperparameter grid and run Hyperopt hyperparameter tuning
    with K-fold cross-validation
    Written by Thejesh Mallidi
    """
    reg = GradientBoostingRegressor(loss=params["loss"],
        learning_rate=params["learning_rate"],
        n_estimators=int(params["n_estimators"]),
        criterion=params["criterion"],
        max_depth=params["max_depth"],
        max_features=params["max_features"],
        random_state=42
    )

    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    validation_loss = cross_validate(
        reg, X_train, y_train,
        scoring="r2",
        cv=cv,
        n_jobs=-1,
        error_score="raise"
    )

    # Note: Hyperopt minimizes the objective, so to maximize R^2, return the negative mean
    return -np.mean(validation_loss["test_score"])


def param_hyperopt(param_grid, max_evals=100):
    """
    Obtain the best parameters from Hyperopt
    Written by Thejesh Mallidi
    """
    trials = Trials()
    params_best = fmin(
        fn=hyperopt_objective,
        space=param_grid,
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials
    )
    
    print("\n\nBest parameters:", params_best)
    return params_best, trials


def gb_reg(trait, fold, n, prefix, plot):
    global X_train
    global y_train

    """ Train Gradient Boosting Regression Model """
    print(trait)

    # Train-test split
    X_train = X.loc[~X.index.isin(test[0])]
    X_test = X.loc[X.index.isin(test[0])]
    y_train = y.loc[~y.index.isin(test[0])]
    y_test = y.loc[y.index.isin(test[0])]

    # Ensure rows are in the same order
    X_train = X_train.loc[y_train.index,:]
    X_test = X_test.loc[y_test.index,:]

    # Hyperparameter tuning
    parameters = {"loss": "absolute_error", # loss function to be optimized
                "learning_rate": hp.uniform("learning_rate", 0.01, 0.5), # learning rate
                "n_estimators": scope.int(hp.quniform("n_estimators", 50, 500, 2),), # number of trees
                "criterion": "friedman_mse", # quality of split
                "max_depth": scope.int(hp.quniform("max_depth", 3, 15, 1)), # tree depth
                "max_features": hp.uniform("max_features", 0.0, 1.0)} # features to consider for split
    
    start = time.time()
    best_params, trials = param_hyperopt(parameters, 100)
    run_time = time.time() - start
    print("Total hyperparameter tuning time:", run_time)
    print("Best parameters: ", best_params)
    print("Trials", trials)

    ################## Training with Cross-Validation ##################
    results_cv = [] # hold performance metrics of cv reps
    results_test = [] # hold performance metrics on test set
    feature_imp = pd.DataFrame(index=X_train.columns)
    preds = {}

    # Training with Cross-validation
    for j in range(0, n): # repeat cv 10 times
        print(f"Running {j+1} of {n}")
        # Build model using the best parameters
        best_model = GradientBoostingRegressor(
            loss="absolute_error",
            learning_rate=best_params["learning_rate"],
            n_estimators=int(best_params["n_estimators"]),
            criterion="friedman_mse",
            max_depth=int(best_params["max_depth"]),
            max_features=best_params["max_features"],
            random_state=j)
        
        cv_pred = cross_val_predict(
            best_model, X_train, y_train, cv=fold, n_jobs=-1) # predictions

        # Performance statistics on validation set
        mse_val = mean_squared_error(y_train, cv_pred)
        rmse_val = np.sqrt(mean_squared_error(y_train, cv_pred))
        evs_val = explained_variance_score(y_train, cv_pred)
        r2_val = r2_score(y_train, cv_pred)
        cor_val = np.corrcoef(np.array(y_train), cv_pred)
        print("Val MSE: %f" % (mse_val))
        print("Val RMSE: %f" % (rmse_val))
        print("Val R-sq: %f" % (r2_val))
        print("Val PCC: %f" % (cor_val[0, 1]))
        result_val = [mse_val, rmse_val, evs_val, r2_val, cor_val[0, 1]]
        results_cv.append(result_val)

        # Evaluate the model on the test set
        best_model.fit(X_train, y_train)
        y_pred = best_model.predict(X_test)

        # Performance on the test set
        mse = mean_squared_error(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        evs = explained_variance_score(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        cor = np.corrcoef(np.array(y_test), y_pred)
        print("Test MSE: %f" % (mse))
        print("Test RMSE: %f" % (rmse))
        print("Test R-sq: %f" % (r2))
        print("Test PCC: %f" % (cor[0, 1]))
        result_test = [mse, rmse, evs, r2, cor[0, 1]]
        results_test.append(result_test)

        # Save the fitted model to a file
        filename = f"{args.save}/{prefix}_model_rep_{j}.pkl"
        pickle.dump(best_model, open(filename, "wb"))

        # Save feature importance scores to file
        feature_imp = pd.concat([feature_imp, pd.Series(best_model.feature_importances_,
            index=best_model.feature_names_in_, name=f"rep_{j}")],
            ignore_index=False, axis=1)
        print(feature_imp.head())

        # Save predicted labels to file
        preds[f"rep_{j}"] = pd.concat([pd.Series(cv_pred, index=X_train.index),
            pd.Series(y_pred, index=X_test.index)], axis=0)

        if plot=="t":
            # Plot linear regression of actual and predicted test values
            sns.regplot(x=y_pred, y=y_test, fit_reg=True, ci=95, seed=123, color="black")
            plt.xlabel("Predicted")
            plt.ylabel("Actual")
            plt.title(f"{trait}")
            plt.savefig(f"{args.save}/{prefix}_lm_test_rep_{j}.pdf", format="pdf")
            plt.close()

    # Write feature importances across reps to file
    feature_imp.to_csv(f"{args.save}/{prefix}_imp.csv")
    
    # Write predictions across reps to file
    pd.DataFrame.from_dict(preds).to_csv(f"{args.save}/{prefix}_preds.csv")
    
    return (results_cv, results_test)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(
        description="Gradient Boosting Regression on SNP and ORF data")
    
    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", help="path to feature table file", required=True)
    req_group.add_argument(
        "-y_name", help="name of label in X file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-save", help="path to save output files", required=True)
    req_group.add_argument(
        "-prefix", help="prefix of output file names", required=True)
    
    # Optional input
    req_group.add_argument(
        "-Y", help="path to label table file", default="")
    req_group.add_argument(
        "-tag", help="description about run to add to results file", default="")
    req_group.add_argument(
        "-fold", help="k number of cross-validation folds", default=5)
    req_group.add_argument(
        "-n", help="number of cross-validation repetitions", default=10)
    req_group.add_argument(
        "-feat", help="file containing features (from X) to include in model", default="all")
    req_group.add_argument(
        "-plot", help="plot feature importances and predictions (t/f)", default="t")
    
    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args() # Read arguments

    # Read in data
    X = dt.fread(args.X)
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)

    if args.Y == "":
        y = X.loc[:, args.y_name]
        X.drop(columns=args.y_name, inplace=True)
    else:
        Y = args.Y
        y = Y.loc[:, args.y_name]
    
    test = pd.read_csv(args.test, header=None)

    # Filter out features not in feat file given - default: keep all
    if args.feat != "all":
        print("Using subset of features from: %s" % args.feat)
        with open(args.feat) as f:
            features = f.read().strip().splitlines()
        X = X.loc[:,features]
        print(f"New dimensions: {X.shape}")

    # Train the model
    start = time.time()
    results_cv, results_test = gb_reg(
        args.y_name, int(args.fold), int(args.n), args.prefix, args.plot)
    run_time = time.time() - start
    print("Training Run Time: %f" % (run_time))

    # Save results to file
    results_cv = pd.DataFrame(
        results_cv, 
        columns=["MSE_val", "RMSE_val", "EVS_val", "R2_val", "PCC_val"])
    results_test = pd.DataFrame(
        results_test, 
        columns=["MSE_test", "RMSE_test", "EVS_test", "R2_test", "PCC_test"])

    # Aggregate results and save to file
    if not os.path.isfile(f"{args.save}/RESULTS_gradboost.txt"):
        out = open(f"{args.save}/RESULTS_gradboost.txt", "a")
        out.write("Date\tRunTime\tTag\tY\tNumInstances\tNumFeatures\
            \tCV_fold\tCV_rep\tMSE_val\tMSE_val_sd\
            \tRMSE_val\tRMSE_val_sd\tEVS_val\tEVS_val_sd\tR2_val\
            \tR2_val_sd\tPCC_val\tPCC_val_sd\tMSE_test\tMSE_test_sd\
            \tRMSE_test\tRMSE_test_sd\tEVS_test\tEVS_test_sd\tR2_test\
            \tR2_test_sd\tPCC_test\tPCC_test_sd")
        out.close()

    out = open(f"{args.save}/RESULTS_gradboost.txt", "a")
    out.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t\
        {run_time}\t{args.tag}\t{args.y_name}\t{X.shape[0]-len(test)}\t\
        {X.shape[1]}\t{int(args.fold)}\t{int(args.n)}\t\
        {np.mean(results_cv.MSE_val)}\t{np.std(results_cv.MSE_val)}\t\
        {np.mean(results_cv.RMSE_val)}\t{np.std(results_cv.RMSE_val)}\t\
        {np.mean(results_cv.EVS_val)}\t{np.std(results_cv.EVS_val)}\t\
        {np.mean(results_cv.R2_val)}\t{np.std(results_cv.R2_val)}\t\
        {np.mean(results_cv.PCC_val)}\t{np.std(results_cv.PCC_val)}\t\
        {np.mean(results_test.MSE_test)}\t{np.std(results_test.MSE_test)}\t\
        {np.mean(results_test.RMSE_test)}\t{np.std(results_test.RMSE_test)}\t\
        {np.mean(results_test.EVS_test)}\t{np.std(results_test.EVS_test)}\t\
        {np.mean(results_test.R2_test)}\t{np.std(results_test.R2_test)}\t\
        {np.mean(results_test.PCC_test)}\t{np.std(results_test.PCC_test)}")
