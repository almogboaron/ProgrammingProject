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
            epsilon = float(arr[2])
            filename1 = arr[3]
            filename2 = arr[4]
        elif len(arr) == 6:
            K = int(arr[1])
            max_iter = int(arr[2])
            epsilon = float(arr[3])
            filename1 = arr[4]
            filename2 = arr[5]
        else:
            raise Exception("Invalid Input!\n")

    except:
        print("Invalid Input!\n")
        sys.exit(1)
    # Asserts
    assert K > 1, "Invalid Input!"
    assert filename1[-4:]==".txt" or filename1[-4:]==".csv" ,"Invalid Input!"
    assert filename2[-4:]==".txt" or filename2[-4:]==".csv" ,"Invalid Input!"
    assert epsilon>0 ,"Invalid Input!"
    assert max_iter>0 ,"Invalid Input!"
    
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
    lst_Indexes = [i for i in range(num_of_rows)]
    lst_choice = [idx]
    distances = np.full(num_of_rows, sys.float_info.max)
    for i in range(1, K):
        for l in range(num_of_rows):
            for j in range(0, i):
                distances[l] = min(distances[l], (np.linalg.norm(np_df[l] - centroids[j]))**2)
        p = np.zeros(num_of_rows)
        for l in range(num_of_rows):
            p[l] = distances[l] / np.sum(distances)
        selected = np.random.choice(lst_Indexes, None, p=p)
        centroids[i] = np_df[selected] # choosing the data frame with max likelihood to be chosen
        lst_choice.append(selected)
    
    #Orginazing Data for C function
    centroids = centroids.tolist()
    np_df = np_df.tolist()

    #Input K,NumOfRows,NumOfCols,max_iter,epsilon,Datalist,CentroidList
    new_Centroids = mykmeanssp.fit(K,num_of_rows,num_of_cols,max_iter,epsilon,np_df,centroids)
    
    #toString
    result = str(lst_choice).strip("[]").replace(" ","")+"\n"
    for m in new_Centroids:
        for i in range(len(m)):
            m[i] = round(m[i],4)
        result += str(m).strip("[]").replace(" ","")+"\n"
    print(result)



if __name__== "__main__":
    main(sys.argv)



