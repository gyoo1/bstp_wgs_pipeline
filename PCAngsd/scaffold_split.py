# grab pandas library
import pandas as pd

# load tsv as dataframe
df = pd.read_csv("scaffold_names_for_split.tsv", sep='\t', header=None)

tot_seq = sum(df[1])

tot_len = len(df[1])

running_count = 0
b_list = []
lengths = []
a_var=0
for x in range(tot_len):
    running_count = running_count + df.iloc[x][1]
    if a_var == 9:
        if x == (tot_len-1):
            lengths.append(running_count)
            b_list.append(x)
    else:
        if running_count > (tot_seq/9):
            b_list.append(x-1)
            lengths.append(running_count - df.iloc[x][1])
            running_count = df.iloc[x][1]
            a_var = a_var + 1

base_str = "region"
counter = 1
sp = 0
for bk in b_list:
    fp = base_str + str(counter) + ".rf"
    if counter == 11:
        break;
    elif counter > 1:
        sp = b_list[counter-2]
        ep = bk
    else:
        ep = bk
    ser = df[sp:ep][0]
    counter = counter + 1
    ser.to_csv(fp,header=False,index=False)