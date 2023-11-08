#!/usr/bin/python3

import numpy as np


def gaussian_grade(B, A, max_score=10, std_dev=1):
    """
    Parameters:
    B (float): The number to be scored.
    A (float): The target number.
    max_score (float): The maximum possible score.
    std_dev (float): The standard deviation of the Gaussian function. A smaller std_dev will
                     make the score drop off more quickly for numbers further from A.

    Returns:
    float: The score of B.
    """
    score = max_score * np.exp(-0.5 * ((B - A) / std_dev) ** 2)
    return score
        
def linear_grade(B, max_score=10, max_value=100):
    """
    Parameters:
    B (float): The number to be scored.
    max_score (float): The maximum possible score.
    max_value (float): The maximum value that B is expected to have. It will get 0 points.

    Returns:
    float: The score of B.
    """
    # 如果B大于最大值，那么得分就是0
    if B >= max_value:
        return 0

    # 否则，得分是线性插值的结果
    score = max_score * (1 - B / max_value)
    return score



