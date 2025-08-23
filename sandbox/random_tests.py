import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt

import matlab.engine as me

eng = me.start_matlab()

# ## construct a bernoulli graph 
# bg = nx.erdos_renyi_graph(100, 0.1)

# config_model = nx.configuration_model([d for n, d in bg.degree()], seed=42)
# degree_dist_bg = bg.degree() 
# degree_dist_config = config_model.degree()



# d_bg = [d for n, d in degree_dist_bg]
# d_config = [d for n, d in degree_dist_config]

# print(d_bg == d_config) ## maybe not the right operation

# # fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# # ax[0].hist(d_bg, bins=10)
# # ax[0].set_title("Degree Distribution (Bernoulli Graph)")

# # ax[1].hist(d_config, bins=10)
# # ax[1].set_title("Degree Distribution (Configuration Model)")

# x = np.arange(1, 10, 0.1)
# y = np.square(x) * np.log(x)
# plt.plot(x, y, label='y = x^2 * log(x)')

# plt.show() 