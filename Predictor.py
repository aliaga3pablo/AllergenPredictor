# -*- coding: utf-8 -*-
"""
Keras sequential model for predicting protein allergenicity

@author: Pablo Aliaga Gaspar
"""
#%%
# Load initial stuff

import pandas as pd
import numpy as np
import tensorflow as tf
import tensorflow_addons as tfa
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate

from imblearn.over_sampling import RandomOverSampler
from tensorflow import keras
from keras.utils import to_categorical
from sklearn.model_selection import train_test_split
from optparse import OptionParser
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve




#%%
def main():
    """
    Description
    -----------
    Optparse function for running the script in BASH
    
    """
    usage = "usage: %prog [options] arg1, arg2, arg3"
    parser = OptionParser(usage)
    parser.add_option("-f", "--file", 
                      action = "store",
                      dest = "filename", 
                      help = "Name of 'data.txt' file with TAB as " 
                      "separator and with format: [First column: Condition], " 
                      "[Second column: Binary condition], [Thir column: Label],"
                      " [Rest: Features]. Use example: filename = 'DF.txt' ")
    parser.add_option("-d", "--dataframe",
                      action = "store_true",
                      dest = "table",
                      default = False,
                      help = "Create a metric data.frame [default: %default]")
    parser.add_option("-g", "--grid",
                      action = "store_true",
                      dest = "grid",
                      default = False,
                      help = "Make a parameter grid for the model with 2 "
                      "hidden layers, trying differents combinations of "
                      "neurons for both layers."
                      "[default: %default]")
    parser.add_option("-j", "--grid2",
                      action = "store_true",
                      dest = "grid2",
                      default = False,
                      help = "Make a parameter grid for the model with 1 "
                      "hidden layer, trying differents combinations of "
                      "neurons and batch size."
                      "[default: %default]")
    parser.add_option("-r", "--roc",
                      action = "store_true",
                      dest = "roc",
                      default = False,
                      help = "Create a ROC curve [default: %default]")
    parser.add_option("-m", "--metrics",
                      action = "store_true",
                      dest = "metrics",
                      default = False,
                      help = "Create 3 curves, each one for 1 metric, "
                      "loss, accuracy and F1 score. "
                      "[default: %default]")
    parser.add_option("-n", "--model",
                      action = "store",
                      dest = "model",
                      default = 2,
                      help = "Select the model to be used, 1 for the model with "
                      "1 hidden layer and 2 for the model with 2 hidden layers. "
                      "[default: %default]")



#%%
# Metrics to be used
METRICS = [
      keras.metrics.BinaryAccuracy(name='accuracy'),
      keras.metrics.Recall(name='recall'),
      keras.metrics.AUC(name='auc'),
      keras.metrics.AUC(name='prc', curve='PR'),
      tfa.metrics.F1Score(num_classes = 2, name="f1")
]

#%%
# Functions definition
def make_model(first, second, metrics=METRICS):
    """
    Description
    -----------
    Keras sequential model definition with a Normalizer layer, 2 Dense hidden
    layers with ReLu activation, 1 Dropout layer and 1 output Dense layer
    with 2 units with Softmax activation.
    The same function will compile the model using Adam as optimizer
    and Categorical Crossentropy as loss metric.
    
    Parameters
    ----------
    first : int
        Number of units for the first hidden layer.
    second : int
        Number of units for the second hidden layer.
    metrics : list, optional
        A list with the metrics that will be measured. The default is METRICS.

    Returns
    -------
    model : keras.engine.sequential.Sequential
        Keras sequential model already defined and compiled.

    """
    model = keras.Sequential([
        normalizer,
        keras.layers.Dense(
            first, activation='relu',
            input_shape=(train_features.shape[-1],)),
        
        keras.layers.Dense(second, activation="relu"),
        keras.layers.Dropout(0.2),
        keras.layers.Dense(2, activation='softmax'),
  ])
  
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=1e-3),
        loss="categorical_crossentropy",
        metrics=metrics)

    return model  


def make_tiny_model(first, metrics=METRICS):
    """
    Description
    -----------
    Keras sequential model definition with a Normalizer layer, a Dense hidden
    layer with ReLu activation, 1 Dropout layer and 1 output Dense layer
    with 2 units with Softmax activation.
    The same function will compile the model using Adam as optimizer
    and Categorical Crossentropy as loss metric.
    
    Parameters
    ----------
    first : int
        Number of units for the hidden layer.
    metrics : list, optional
        A list with the metrics that will be measured. The default is METRICS.

    Returns
    -------
    model : keras.engine.sequential.Sequential
        Keras sequential model already defined and compiled.

    """
    model = keras.Sequential([
        normalizer,
        keras.layers.Dense(
            first, activation='relu',
            input_shape=(train_features.shape[-1],)),
        
        keras.layers.Dropout(0.2),
        keras.layers.Dense(2, activation='softmax'),
  ])
  
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=1e-3),
        loss="categorical_crossentropy",
        metrics=metrics)

    return model  


def parameter_grid(values_list):
    """
    Desciption
    ----------
    Parameter grid to use with the model with 2 hidden layers, this
    function will determine the best combination of units to maximize
    F1 score for test set.
    
    Parameters
    ----------
    values_list : list
        Range of values to be tested as hidden neurons.

    Returns
    -------
    f1_df : pandas.DataFrame
        DataFrame that contains F1 score for test set and the number of
        units for each layer.

    """
    f1_list = []
    range_to_test = values_list
    
    for i in range_to_test:
        for j in range_to_test:
            print("Running... Fist Layer ",j, " and Second Layer: ",i)
            model = make_model(j,i)
            
            baseline_history = model.fit(
                train_features,
                train_labels,
                batch_size=BATCH_SIZE,
                epochs=EPOCHS,
                callbacks=[early_stopping],
                validation_data=(val_features, val_labels))
            
            results = model.evaluate(test_features, test_labels)
            f1_list.append(results[6].mean())
    f1_matrix = np.array(f1_list).reshape(int(np.sqrt(len(f1_list))),int(np.sqrt(len(f1_list))))
    f1_df = pd.DataFrame(f1_matrix)
    f1_df = f1_df.set_axis(range_to_test, axis=1)
    f1_df = f1_df.set_axis(range_to_test, axis=0)
    
    print("First layer for best F1: " + str(f1_df.max(axis=1).idxmax()))
    print("Second layer for best F1: " + str(f1_df.max(axis=0).idxmax()))
    print("Best F1 score: " + str(
        f1_df.loc[f1_df.max(axis=1).idxmax(), f1_df.max(axis=0).idxmax()]))


    ax = sns.heatmap(f1_df, annot=False)
    ax.set(xlabel = "First layer neurons", ylabel = "Second layer neurons")
    plt.savefig('heatmap.png')
    plt.show()
    
    return f1_df
    

def parameter_grid_2(values_list, batch_list):
    """
    Description
    -----------
    Parameter grid to use with the model with 1 hidden layer, this
    function will determine the best combination of units and batch size
    to maximize F1 score for test set.

    Parameters
    ----------
    values_list : list
        Range of values to be tested as hidden neurons.
    batch_list : list
        Range of values to be tested as batch size.

    Returns
    -------
    f1_df : pandas.DataFrame
        DataFrame that contains F1 score for test set and the number of
        units the layer and the batch size.

    """
    f1_list = []
    range_to_test = values_list
    batch_to_test = batch_list

    for i in range_to_test:
        for j in batch_to_test:
            print("Running... Layer: ",i, " and Batch size: ",j)
            model = make_tiny_model(i)
            
            baseline_history = model.fit(
                train_features,
                train_labels,
                batch_size=j,
                epochs=EPOCHS,
                callbacks=[early_stopping],
                validation_data=(val_features, val_labels))
            
            results = model.evaluate(test_features, test_labels)  
            f1_list.append(results[6].mean())
    f1_matrix = np.array(f1_list).reshape(len(range_to_test), len(batch_to_test))
    f1_df = pd.DataFrame(f1_matrix)
    f1_df = f1_df.set_axis(range_to_test, axis=0)
    f1_df = f1_df.set_axis(batch_to_test, axis=1)

    print("Batch size for best F1: " + str(f1_df.max(axis=1).idxmax()))
    print("Layer units for best F1: " + str(f1_df.max(axis=0).idxmax()))
    print("Best F1 score: " + str(
        f1_df.loc[f1_df.max(axis=1).idxmax(), f1_df.max(axis=0).idxmax()]))


    ax = sns.heatmap(f1_df, annot=False)
    ax.set(xlabel = "Batch size", ylabel = "Hidden Layer neurons")
    plt.savefig('heatmap_tiny.png')
    plt.show()  
    
    return f1_df


def make_roc():
    """
    Description
    -----------
    Function to make a ROC curve using the test set.

    Returns
    -------
    Print ROC curve

    """
    predicted_labels = model.predict(test_features)[:,1]

    auc = roc_auc_score(test_labels[:,1], predicted_labels)
    pred_fpr, pred_trp, _ = roc_curve(test_labels[:,1], predicted_labels)

    plt.plot(pred_fpr, pred_trp, label = f"AUC = {auc.round(3)}")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend()
    plt.savefig('roc_curve.png')
    

def plot_metrics():
    """
    Description
    -----------
    Function to print metrics for training and validation sets during
    each training epoch.
    
    Returns
    -------
    Print Loss, Accuracy and F1 curves
    """
    # Loss
    ax_loss = sns.lineplot(data=history_df[["loss", "val_loss"]])
    ax_loss.set(xlabel = "Epochs", ylabel = "Loss")
    plt.savefig('loss_curve.png')
    plt.show()
    # Accuracy
    ax_accuracy = sns.lineplot(data=history_df[["accuracy", "val_accuracy"]])
    ax_accuracy.set(xlabel = "Epochs", ylabel = "Accuracy")
    plt.savefig('accuracy_curve.png')
    plt.show()
    # F1 score
    ax_f1 = sns.lineplot(data=history_df[["f1", "val_f1"]])
    ax_f1.set(xlabel = "Epochs", ylabel = "F1 score")
    plt.savefig('f1_curve.png')
    plt.show()
    

#%%
# Run optparse function
parser, options = main()

# Dataset preprocessing
# Load dataset
DF = pd.read_csv(options.filename, sep = "\t")

# Transform columns with big scale differences to logarithmic
DF["log_molecular_weight"] = np.log(DF.pop("molecular_weight"))
DF["log_pepLength"] = np.log(DF.pop("pepLength"))

# Split dataset into training, validation and testing
train_df, test_df = train_test_split(DF, test_size = 0.2, random_state = 123)
train_df, val_df = train_test_split(train_df, test_size = 0.2, random_state = 123)

# Labels separation
train_labels = to_categorical(np.array(train_df.pop("binCondition")))
val_labels = to_categorical(np.array(val_df.pop("binCondition")))
test_labels = to_categorical(np.array(test_df.pop("binCondition")))

# Remove unnecesary columns
train_features = np.array(train_df.drop(["condition", "protNames"], axis=1))
val_features = np.array(val_df.drop(["condition", "protNames"], axis=1))
test_features = np.array(test_df.drop(["condition", "protNames"], axis=1))

# Random oversampling for solving umbalancement
oversample = RandomOverSampler(sampling_strategy='minority')

train_features, train_labels = oversample.fit_resample(train_features, train_labels)
train_labels = to_categorical(train_labels)


#%%
# Normalization layer definition
normalizer = keras.layers.Normalization()
normalizer.adapt(train_features)


#%%
# Parameters that will be used when building de model
EPOCHS = 1000
BATCH_SIZE = 1024

early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss', 
    verbose=1,
    patience=10,
    mode='min',
    restore_best_weights=True)


#%%
# Parameter grid
if options.grid is True:
    range_to_test = range(5,205,5)
    f1_df = parameter_grid(range_to_test)
    with open("f1_DF.txt", "w") as f:
        f.write(tabulate(f1_df, headers = f1_df.columns.values,
                         tablefmt="tsv"))

#%%
# Parameter grid for tiny model
if options.grid2 is True:
    range_to_test = range(5,400,5)
    batch_to_test = range(30,2000,10)
    f1_tiny_df = parameter_grid_2(range_to_test, batch_to_test)
    with open("f1_tiny_DF.txt", "w") as f:
        f.write(tabulate(f1_tiny_df, headers = f1_tiny_df.columns.values,
                         tablefmt="tsv"))

#%%
# Model creation 
if options.model == 2:
    model = make_model(160, 190)
elif options.model == 1:
    model = make_tiny_model(60, 380)

model.summary()


history = model.fit(
    train_features,
    train_labels,
    batch_size=BATCH_SIZE,
    epochs=EPOCHS,
    callbacks=[early_stopping],
    validation_data=(val_features, val_labels))


# Obtain mean for F1 score interval
history_df = pd.DataFrame.from_dict(history.history)

f1_train = []
for element in history_df["f1"]:
    f1_train.append(element.mean())
    
f1_val = []
for element in history_df["val_f1"]:
    f1_val.append(element.mean())

history_df["f1"] = f1_train
history_df["val_f1"] = f1_val


#%%
# Plot metrics (Loss, Accuracy and F1 score for training and validation sets)
if options.metrics is not False:
    plot_metrics()

# Receiver operating characteristics (ROC) curve and area under the curve (AUC)
if options.roc is not False:
    make_roc()


# Model evaluation with test data
test_scores = model.evaluate(test_features, test_labels)
results_df = pd.DataFrame(test_scores)
results_df.columns = ["Metric score"]
results_df.index = model.metrics_names
results_df.iloc[5,0] = results_df.iloc[5,0].mean()

if options.table is not False:
    with open("resultsDF.txt", "w") as f:
        f.write(tabulate(results_df, headers = results_df.columns.values,
                         tablefmt="tsv"))
        
#%%
# Table to identify new possible allergens
predictions = model.predict(test_features)[:,1]
predictions_df = test_df.iloc[:,0:2].assign(Prediction=predictions)
predictions_df = predictions_df.sort_values(by="Prediction", ascending=False)
predictions_df = predictions_df[predictions_df["condition"] == "Non allergen"]


with open("predictions_DF.txt", "w") as f:
    f.write(tabulate(predictions_df, headers = predictions_df.columns.values,
                     tablefmt="tsv"))



