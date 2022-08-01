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
def kmeans( matrix ):
    # kmeans++ implementation
    matrix = np.matrix(matrix)
    centroids = np.zeros(shape=(matrix.shape[1], matrix.shape[1]))
    idx = np.random.choice(matrix.shape[0])  # choosing random row index
    centroids[0] = matrix[idx]
    lst_Indexes = [i for i in range(matrix.shape[0])]
    lst_choice = [idx]
    distances = np.full(matrix.shape[0], sys.float_info.max)
    for i in range(1, matrix.shape[1]):
        for l in range(matrix.shape[0]):
            for j in range(0, i):
                distances[l] = min(distances[l], (np.linalg.norm(matrix[l] - centroids[j])) ** 2)
        p = np.zeros(matrix.shape[0])
        for l in range(matrix.shape[0]):
            p[l] = distances[l] / np.sum(distances)
        selected = np.random.choice(lst_Indexes, None, p=p)
        centroids[i] = matrix[selected]  # choosing the data frame with max likelihood to be chosen
        lst_choice.append(selected)

    # Orginazing Data for C function
    centroids = centroids.tolist()
    matrix = matrix.tolist()

    # Input K,NumOfRows,NumOfCols,max_iter,epsilon,Datalist,CentroidList
    new_Centroids = spkmeansmodule.kmeans_fit(matrix.shape[1], matrix.shape[0], matrix.shape[1], 300, 0, matrix, centroids)

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

    # Numpy Random Seed :
    np.random.seed(0)

    # Cases Functions
    if(goal.val == 3):
        newCentroids = spkmeansmodule.spk()
        kmeans(newCentroids)

    elif(goal.val == 2):
        spkmeansmodule.wam()

    elif(goal.val == 3):
        spkmeansmodule.ddg()

    elif(goal.val == 4):
        spkmeansmodule.lnorm()

    elif(goal.val == 5):
        spkmeansmodule.jacobi()


if __name__== "__main__":
    main(sys.argv)