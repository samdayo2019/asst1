import subprocess
import pandas as pd
import os


NUM_RUNS = 10

def main():

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    total_assignments = 0
    total_centroids = 0
    total_cost = 0
    total_time = 0

    for i in range(NUM_RUNS):
        subprocess.run(["./kmeans"], check=True)

        file_name = 'Average_Profiled_Elapsed_Times.csv'
        file_name2 = 'Total_kmeans_time.csv'

        if(os.path.exists(file_name) and os.path.exists(file_name2)):
            df_processing = pd.read_csv(file_name)
            df_processing_time = pd.read_csv(file_name2)

            assignments_val = df_processing['Assignment'].iloc[0]
            centroids_val = df_processing['Centroids'].iloc[0]
            cost_val = df_processing['Cost'].iloc[0]
            elapsed_time = df_processing_time['Run time for Kmeans'].iloc[0]

            total_assignments += assignments_val
            total_centroids += centroids_val
            total_cost += cost_val
            total_time += elapsed_time

        else: 
            print("CSV files not found after iterations {}".format(i+1))
            continue


    average_assignments = total_assignments/NUM_RUNS
    average_centroids = total_centroids/NUM_RUNS
    average_cost = total_cost/NUM_RUNS
    average_time = total_time/NUM_RUNS

    print("Average Assignments: {}".format(average_assignments))
    print("Average Centroids: {}".format(average_centroids))
    print("Average Cost: {}".format(average_cost))
    print("Average Time: {}".format(average_time))



if __name__== "__main__":
    main()  