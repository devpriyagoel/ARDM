########################
# Author: Dev Priya Goel
########################


import numpy as np
import math

def calc_numerical_gradient(func, x, delta_x):
    """Function for computing gradient numerically."""
    val_at_x = func(x)
    val_at_next = func(x + delta_x)
    return (val_at_next - val_at_x) / delta_x


def ardm(func, L, dimension, init_x=None, numerical_gradient=True, delta_x=0.0005, gradient_func=None,epsilon=None):
    
    assert delta_x > 0, "Step must be positive."

    if (init_x is None):
        x = np.zeros(dimension)  # todo проверка что функция определена в 0.
    else:
        x = init_x

    if (epsilon is None):
        epsilon = 0.05

    alpha = 0.05 / (2 * L)
    prev_x = np.zeros(dimension)
    v = np.zeros(dimension)
    temp_ans = np.zeros(dimension)
    if numerical_gradient:
        gradient = calc_numerical_gradient(func, x, delta_x)
    else:
        gradient = gradient_func(x)
    gc = np.linalg.norm(gradient)
    gp = 0.0
    iterations =0
    max_iterations = 10000
    step=0
    while iterations<max_iterations:
        iterations+=1
        if step==0:
            v = x - alpha*gradient
            prev_x = x
            if numerical_gradient:
                temp_ans = calc_numerical_gradient(func, v, delta_x)
            else:
                temp_ans = gradient_func(v)
            x = v-alpha*temp_ans
        else:
            beta = gc/gp
            v = x + beta*(x-prev_x) - alpha*(1+beta)*gradient
            prev_x = x
            if numerical_gradient:
                temp_ans = calc_numerical_gradient(func, v, delta_x)
            else:
                temp_ans = gradient_func(v)
            x = v-alpha*temp_ans
        gp = gc
        if numerical_gradient:
            gradient = calc_numerical_gradient(func, x, delta_x)
        else:
            gradient = gradient_func(x)
        gc = np.linalg.norm(gradient)
        if gc>gp&step!=0:
            uc = up
            k=-1
        if gc<epsilon&step!=0:
            break

    return x
