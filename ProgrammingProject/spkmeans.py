import numpy as np
import pandas as pd
import random
import sys
import enum
import spkmeansmodule

# Enum Class For Goals
class Goal(enum.Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5

# Kmeans++ Func for (include kmeans through C)
def kmeans(matrix):
    # kmeans++ implementation
    matrix = np.matrix(matrix)
    K = matrix.shape[1]
    centroids = np.zeros(shape=(K, matrix.shape[1]))
    idx = np.random.choice(matrix.shape[0])  # choosing random row index
    centroids[0] = matrix[idx]
    lst_Indexes = [i for i in range(matrix.shape[0])]
    lst_choice = [idx]
    distances = np.full(matrix.shape[0], sys.float_info.max)
    for i in range(1, K):
        for l in range(matrix.shape[0]):
            for j in range(0, i):
                distances[l] = min(distances[l], (np.linalg.norm(matrix[l] - centroids[j])) ** 2)
        p = np.zeros(matrix.shape[0])
        for l in range(matrix.shape[0]):
            p[l] = distances[l] / np.sum(distances)
        selected = np.random.choice(lst_Indexes, None, p=p)
        centroids[i] = matrix[selected]  # choosing the data frame with max likelihood to be chosen
        lst_choice.append(selected)

    # Organizing Data for C function
    centroids = centroids.tolist()

    # CentroidList
    new_Centroids = spkmeansmodule.fit(centroids)

    # toString
    result = str(lst_choice).strip("[]").replace(" ", "") + "\n"
    for m in new_Centroids:
        for i in range(len(m)):
            m[i] = float("{:.4f}".format(m[i]))
        result += str(m).strip("[]").replace(" ", "") + "\n"
    print(result)

# Main Func
def main(arr):

    # reading user CMD
    try:
        if len(arr) == 4:
            K = int(arr[1])
            goal = Goal[arr[2].lower()]
            file_name = arr[3]
        else:
            raise Exception("Invalid Input!\n")
    except:
        print("Invalid Input!\n")
        sys.exit(1)

    # Asserts :
    assert file_name[-4:] == ".txt" or file_name[-4:] == ".csv", "Invalid Input!"
    assert K >= 0, "Invalid Input!"

    # Numpy Random Seed :
    np.random.seed(0)

    # Cases Functions
    if goal.value == 1:
        U = spkmeansmodule.spkC(file_name, K)
        kmeans(U)

    elif goal.value == 2:
        spkmeansmodule.wamC(file_name)

    elif goal.value == 3:
        spkmeansmodule.ddgC(file_name)

    elif goal.value == 4:
        spkmeansmodule.lnormC(file_name)

    elif goal.value == 5:
        spkmeansmodule.jacobiC(file_name)


if __name__ == "__main__":
    main(sys.argv)