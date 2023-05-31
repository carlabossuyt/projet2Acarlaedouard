import pandas as pd
import matplotlib.pyplot as plt


data=[[0,-2,-4,-6,-8],
      [-2,0,0,0,0],
      [-4,0,0,0,0],
      [-6,0,0,0,0],
      [-8,0,0,0,0]]
column_labels=["","A", "T", "C", "G"]
df=pd.DataFrame(data,columns=column_labels)
plt.axis('tight')
plt.axis('off')
plt.table(cellText=df.values,colLabels=df.columns,rowLabels=["","G","C","T","A"],loc="center")

plt.show()
