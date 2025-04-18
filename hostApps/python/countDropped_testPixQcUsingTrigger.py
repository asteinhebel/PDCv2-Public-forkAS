import pandas as pd
import sys
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1], sep=';')

#Get trigger rows
expectedTrigs = 5
df_trigs = df[[f'SPAD_TCR{i}' for i in range(4)]]
df_droppedTrigs = df_trigs[df_trigs.isin([4]).any(axis=1)]

print(f"Number of QCs that did not record {expectedTrigs} triggers")
print((df_droppedTrigs != expectedTrigs).sum())

print(f"Index locations of QCs that not record {expectedTrigs} triggers")
print(df_droppedTrigs.head(-1))

#plot
plt.figure()
s=[100,50,30,10]
for i in range(4):
    alpha=0.6 if i>1 else 1
    plt.scatter(df_droppedTrigs.index,df_droppedTrigs[f'SPAD_TCR{i}'], label=i, s=s[i])
plt.legend(loc='best')
plt.show()
plt.close()




