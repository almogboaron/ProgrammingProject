import numpy as np
import pandas as pd
import random
import sys
import mykmeanssp



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

    num_of_cols = df.shape[1]
    num_of_rows = df.shape[0]

    # kmeans++ implementation
    centroids = np.zeros(shape=(K,num_of_cols))
    np_df = df.to_numpy() # working with our dataframe in numpy format
    np.random.seed(0)
    idx = np.random.choice(num_of_rows) # choosing random row index
    centroids[0] = np_df[idx]
    lst_Indexes = [idx]
    for i in range(1, K):
        distances = np.full(num_of_rows, sys.float_info.max)
        distances[0] = np.linalg.norm(np_df[0] - centroids[0])
        for l in range(num_of_rows):
            for j in range(1, i):
                distances[l] = min(distances[l], np.linalg.norm(np_df[l] - centroids[j]))
        p = np.zeros(num_of_rows)
        for l in range(num_of_rows):
            p[l] = distances[l] / np.sum(distances)
        centroids[i] = np_df[np.argmax(p)] # choosing the data frame with max likelihood to be chosen
        lst_Indexes.append[np.argmax(p)]
    #Orginazing Data for C function
    centroids = np.tolist(centroids)
    np_df = np.tolist(np_df)

    #Input K,NumOfRows,NumOfCols,max_iter,epsilon,Datalist,CentroidList
    new_Centroids = mykmeanssp.fit(K,num_of_rows,num_of_cols,max_iter,epsilon,np_df,centroids)
    print(lst_Indexes)
    for centroid in new_Centroids:
        print(centroid)


if __name__== "__main__":
    main(sys.argv)



