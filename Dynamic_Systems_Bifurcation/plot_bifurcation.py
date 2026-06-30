import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("euler_b.csv")

plt.figure(figsize=(12, 8))

plt.scatter(
    data["parameter"],
    data[" peak"],
    s=1
)

plt.xlabel("b")
plt.ylabel("Peak x")
plt.title("Bifurcation diagram")

plt.show()