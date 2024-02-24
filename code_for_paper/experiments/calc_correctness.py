import numpy as np

# The parameters should be the pandas dataframe
def calc_correctness(PC1_df, approx_df):
    PC1_np = PC1_df.values # Turn into numpy format
    PC1_np = PC1_np.flatten() # Turn into 1D vector
    PC1_np = PC1_np[PC1_np != 0] # Remove 0

    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector
    approx_np = approx_np[approx_np != 0] # Remove 0

    del PC1_df, approx_df

    if len(PC1_np) != len(approx_np): 
        print("PC1 and approx has a different number of elements")
        return
    
    valid_entry_num = len(PC1_np)

    if np.corrcoef(PC1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    PC1_pos_np = PC1_np > 0
    approx_pos_np = approx_np > 0
    PC1_pos_VS_approx_pos_np = PC1_pos_np == approx_pos_np 
    
    correct_num = list(PC1_pos_VS_approx_pos_np).count(True)
    correct_rate = correct_num / valid_entry_num

    return {
        "valid_entry_num": valid_entry_num,
        "correct_num": correct_num,
        "correct_rate": correct_rate,
    }