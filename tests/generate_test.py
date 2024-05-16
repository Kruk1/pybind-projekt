import scikit_build_example as m
import numpy as np

n = np.arange(2000)
sin = m.generate_signal(100, 48000, n, 1, 0)

m.visualization(sin)

cos = m.generate_signal(100, 48000, n, 2, 0)

m.visualization(cos)

pilo = m.generate_signal(100, 48000, n, 3, 0)

m.visualization(pilo)

rectangle = m.generate_signal(100, 48000, n, 4, 0)

m.visualization(rectangle)