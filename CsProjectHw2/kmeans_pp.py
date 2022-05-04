import numpy as np
import pandas as pd
import random
import sys



def main(arr):

    # reading user CMD
    try:
        if len(arr) == 5:
            K = int(arr[1])
            max_iter = 300
            epsilon = arr[2]
            filename1 = arr[3]
            filename2 = arr[4]
        elif len(arr) == 6:
            K = int(arr[1])
            max_iter = int(arr[2])
            epsilon = arr[3]
            filename1 = arr[4]
            filename2 = arr[5]
        else:
            raise Exception("Invalid Input!\n")

    except:
        print("Invalid Input!\n")
        sys.exit(1)

    assert K > 1, "Invalid Input!"

    # inner join of 2 files
    table1 = pd.read_csv(filename1,header=None)
    table2 = pd.read_csv(filename2,header=None)
    df = pd.merge(table1,table2,on=0,how="inner",sort=True)
    df.drop(columns=df.columns[0],axis=1,inplace=True)

    # kmeans++ implementation
    centroids = np.zeros(shape=(K,len(df.columns)))
    np.random.seed(0)
    centroids[0] = np.random.choice(df.rows)
    for i in range(1,K):
        distances = np.linalg.norm(df.to_numpy() - centroids[0])
        for j in range(1, i):
            distances = np.minimum(distances, np.linalg.norm(df.to_numpy() - centroids[j]))










if __name__== "__main__":
    main(sys.argv)



