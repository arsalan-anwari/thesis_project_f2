# Dynamic Training Routine 
This documents describes the different settings that can be tweak in the simulation code of the `dynamic-training-routine.R` script.

## Settings
There a number of global variables with suffix `settings._` which influence the behaviour of the training routine but the only one that should be modified is `settings.loss_function.default_power`. This is a scalar value between `0.0` and `1.0` with a default of `0.68732`. The lower this value the more drastic loss values calculates by the loss function will be.

## Weight functions
These functions with suffix `weight_func._` are responsible for the translation of loss values to a correct hyper parameter. They have the following prototype:

`weight_func._(old_weight, loss_value, weight_func_opts = FALSE) -> new_weight`
- `{old:new}_weight := any`
- `loss_value := float`
- `weight_func_opts := float or list(params...)`

It possible to define your own weight function following the prototype listed above, but some weight functions are already defined these being:

| Description | Input/Output |
|-|-|
| `one_base_sub_halve (...)` <li> Adds one by default. <li> Increments by `loss_value`. <li> Decrements by halve of `loss_value` by default <li> Increments by `loss_value` if weights becomes negative      | <li> `IN := {0.0 ... 1.0}` <li> `OUT := {1.0 ... n}` |
| `no_base_add_exp_sub_halve (..., exp) ` <li> Adds nothing by default. <li> Increments by `loss_value^exp`. <li> Decrement method same as `one_base_sub_halve` | <li> `IN := {0.0 ... 1.0}` <li> `OUT := {0.0 ... n}` |
| `rand_base_add_sub_inc` <br> `(..., list(n_weights: int , alpha: float)) ` <li> Adds `runif(n_weights, 0+alpha, 1-alpha)` by default<li> Increments by `loss_value`. <li> Decrement method by `loss_value + alpha` | <li> `IN := {0.0 ... 1.0}` <li> `OUT := {0.0 ... n}[n_weights]` |
| `variogram_base_add_exp` <br> `(..., list(sill_power: int , range_power: int, exp: float)) ` <li> multiplies result of `no_base_add_exp_sub_halve(weight)` with `(psill, range)` | <li> `IN := list(weight, psill, range)` <li> `OUT := list(weight, psill, range)` |

## RMSE functions

The dynamic training model requires a valid function which calculates the RMSE of a given model using a cross validation method and has the following prototype:

`rmse_func._(weights, model_opts = FALSE) -> {0.0 ... n}`
- `weights := any`
- `model_opts := any`

It possible to define your own RMSE function following the prototype listed above, but RSME functions using SLOO-CV are already defined these being:
| RSME Function | Params
|-|-|
| `rmse_func.trend_surface` | <li> `case_weight := {1.0 ... n}` <li> `formula := string` | 
| `rmse_func.idw` | <li> `idp_power := {1.0 ... n}` <li> `neighbors := int` | 
| `rmse_func.od_kriging` | <li> `vgm_weights := list(weight: float, psill: float, range: float)` <li> `neighbors := int` | 
| `rmse_func.u_kriging` | <li> `vgm_weights := list(weight: float, psill: float, range: float)` <li> `model_opts := list(c(int, string), ...)` | 


## Calculating model parameters

The dynamic training routine can be started calling the function `dynamic_model.train_opts(...) -> data.frame()` with the following settings:
```SQL
  df_res <- dynamic_model.train_opts(
    rmse_func: function(any, ?any) -> float,
    model_opts: ?any,
    weight_func: function(any, float, ?any) -> any,
    weight_func_opts: ?any,
    max_iterations: ?int = 20,
    loss_value_threshold: ?float,
    loss_function_power: ?float = settings.settings.loss_function.default_power 
  )

```

- `max_interation`: maximum number of iterations to run each variation of the model
- `loss_value_threshold`: threshold value where the routine is stopped when calculated `loss_value < loss_value_threshold`. This will stop the routine prematurely even if `i < max_iterations`

## Output format

The output format of `df_res` will be a `data.frame()` with the following columns.

```SQL
model_opt: int
train_iteration: int
old_rsme: float
new_rsme: float
loss_value: float
weights: char
```

When exporting `df_res` to a `.rds` file all columns will become strings by default, so some type conversion needs to be peformed before it can be used in another script or notebook. 

The column `weights` contains the user defined hyperparameters in the format of `model_opts` decoded as string with format `"(model_opts[1], model_opts[n], ...)"`. 
