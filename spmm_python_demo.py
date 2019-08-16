import numpy as np

val = np.array([1/3, 1/2, 1/2, 1/3, 1/2, 1, 1/2, 1/3, 1/3, 1/3, 1/3, 1/2, 1/3, 1/2, 1/3, 1, 1/3])
col_idx = np.array([6, 0, 2, 3, 0, 1, 2, 3, 6, 3, 4, 7, 4, 7, 4, 5, 6])
row_ptr = np.array([0, 1, 4, 5, 6, 9, 12, 14])
x0 = (1/8)*np.ones(8)
C = np.zeros(len(row_ptr))

for i in range(len(row_ptr)):
    for j in range(row_ptr[i], row_ptr[i+1] - 1):
        C[i] += val[j]*x0[int(col_idx[j])]

"""for j in range(row_ptr[-1], len(val)):
    C[-1] += val[j]*x0[int(col_idx[j])]"""

print(C)

# eureka
