import sys
import numpy as np
import pandas as pd
import mykmeanssp

try:
    k = int(sys.argv[1])
    if k <= 0:
        print("Error in K argument!")
        exit()
except ValueError:
    print("Error in K argument!")
    exit()
except IndexError:
    print("No arguments provided!")
    exit()

max_iter = 300
if len(sys.argv) == 5:
    try:
        max_iter = int(sys.argv[2])
        if max_iter < 0:
            print("Error in max_iter argument!")
            exit()
    except ValueError:
        print("Error in max_iter argument!")
        exit()

    input1_pd = pd.read_csv(sys.argv[3], header=None)
    input2_pd = pd.read_csv(sys.argv[4], header=None)
else:
    input1_pd = pd.read_csv(sys.argv[2], header=None)
    input2_pd = pd.read_csv(sys.argv[3], header=None)

df = pd.merge(input1_pd, input2_pd, on=0)
df = df.sort_values(df.columns[0])
df = df.reset_index(drop=True)
df = df.drop(df.columns[0], axis=1)

dim = len(df.columns)

df['cluster'] = -1.0
points_to_cluster_init = df.to_numpy().flatten().tolist()
del df['cluster']

np_df = df.to_numpy()
num_of_points = len(np_df)

df['d'] = np.nan
df['p'] = np.nan
np_df_ext = df.to_numpy()

if k >= len(df):
    print("Error in K argument!")
    exit()

np.random.seed(0)
rand_start = np.random.randint(0, (len(df) - 1))
centroids_df_indices = [rand_start]
z = 1

while z < k:
    for i in range(num_of_points):
        ind = centroids_df_indices[z - 1]
        d_cand = (np.linalg.norm(np_df[i] - np_df[ind])) ** 2
        if z == 1:
            np_df_ext[i][dim] = d_cand
        else:
            if d_cand < np_df_ext[i][dim]:
                np_df_ext[i][dim] = d_cand
    d_sum = np_df_ext.sum(axis=0)[dim]
    for i in range(num_of_points):
        np_df_ext[i][dim + 1] = np_df_ext[i][dim] / d_sum

    n_array = np.arange(num_of_points)
    centroids_df_indices.append(np.random.choice(n_array, p=np_df_ext[:, dim + 1]))
    z += 1

centroids_lists = [np_df[ind].tolist() for ind in centroids_df_indices]  # need to output this to file
centroids_init = [num for sublist in centroids_lists for num in sublist]

num_of_lines = len(df)
to_print = ','.join([str(num) for num in centroids_df_indices])
print(to_print)

centroids = mykmeanssp.fit(k, max_iter, num_of_lines, dim, centroids_init, points_to_cluster_init,
                           len(centroids_init),
                           len(points_to_cluster_init))

centroids = np.asarray(centroids).round(4)
centroids = np.array_split(centroids, k)
centroids = [centroids[i].tolist() for i in range(len(centroids))]

for cent in centroids:
    print(','.join(str(val) for val in cent))
