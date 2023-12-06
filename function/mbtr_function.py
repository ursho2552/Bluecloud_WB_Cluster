def mbtr_fit(path, hp_id, loss_type: str='mse', n_boosts: int=1000, n_q: int=10, learning_rate: float=0.1, 
early_stopping_rounds: int= 3,
min_leaf: int=20, lambda_weights: float=0.0001, lambda_leaves: float=0.0001):

  """
  Fit a MBTR model under mean square error loss function. The inputs corresponds
  to the MBT().fit() function of the MBTR module except:
    
  :input x_tr: training set generated by the 05_folds.R script
  :input y_tr: training set generated by the 05_folds.R script
  
  (optionnal)
  :input x_val: external validation set generated by the 05_folds.R script
  :input y_val: external validation set generated by the 05_folds.R script
  
  This function returns a list containing:
  [[1]] the MBT model object
  [[2]] a list of the loss computed on the validation datasets for each boosting round
  Results are saved in the cv_hp_m object and further concatenated in a list in R.
  
  See the MBT().fit() and https://mbtr.readthedocs.io/en/latest/mbtr.html#module-mbtr.mbtr
  for more details on parameters.
  """

  import numpy as np
  import mbtr.utils as ut
  from mbtr.mbtr import MBT
  from mbtr.utils import set_figure
  from pandas import read_feather
  import pickle
  
  # Load X and Y from R feather file
  X_tr = read_feather(path+"X_tr.feather").to_numpy()
  Y_tr = read_feather(path+"Y_tr.feather").to_numpy()
  
  # Load X_val and Y_val from R environment
  X_val = read_feather(path+"X_val.feather").to_numpy()
  Y_val = read_feather(path+"Y_val.feather").to_numpy()
  
  # Fit
  m = MBT(loss_type = loss_type,
          n_boosts=n_boosts, 
          n_q=n_q,
          learning_rate=learning_rate, 
          min_leaf=min_leaf,
          lambda_weights=lambda_weights, 
          lambda_leaves=lambda_leaves,
          verbose = 0,
          early_stopping_rounds=early_stopping_rounds).fit(X_tr, Y_tr, X_val, Y_val, do_plot=False)
          
  # Save model        
  with open(path+hp_id+'_m', 'wb') as f:
    pickle.dump(m, f, pickle.HIGHEST_PROTOCOL)
          
  return m

def mbtr_predict(model, X_pred, n_boosts):

  """
  Predict the y_hat target dataframe according to a x_pred feature matrix and
  a previously fitted MBTR model.
  
  :param model: a MBT model object
  :param X_pred: a Nfeature * Nobs matrix on which the model will predict the
                 target values y_hat
                 
  See the MBT().fit() and https://mbtr.readthedocs.io/en/latest/mbtr.html#module-mbtr.mbtr
  for more details on parameters.
  """

  import numpy as np
  import mbtr.utils as ut
  from mbtr.mbtr import MBT
  from mbtr.utils import set_figure
  
  # Load X_pred from R environment
  X_pred = X_pred.to_numpy()
  
  # Predictions
  y_hat = model.predict(x = X_pred, n = n_boosts)
  
  return y_hat
  
  
