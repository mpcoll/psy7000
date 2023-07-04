# TEST
Créé par Michel-Pierre Coll - 2022, michel-pierre.coll@psy.ulaval.ca


## 0. Instructions - calepin interactifs

Les simulations, la visualisation et l'interaction sont d'excellentes façon de comprendre certains concepts statistiques (Moreau, 2015). Au cours de la session nous utiliserons régulièrement des simulations et visualisations pour illuster certains concepts. 

Ces simulations sont effectuées en utilisant des documents comme celui-ci ("calepins" ou *notebooks*) qui contiennent des explications et des "cellules" qui permettent d'exécuter le langage de programmation *Python*. Dans la majorité des cas, le code sera caché. Il vous sera parfois demandé de changer quelques paramètres dans le code pour évaluer leur effet sur les résultats à l'aide de boutons. **LA COMPRÉHENSION DU CODE DANS LES CELLULES N'EST PAS NÉCESSAIRE POUR LE COURS.** Si vous aimeriez comprendre comment les simulations sont effectuées, vous pouvez cliquer sur "Show code".


Pour débuter, exécutez le calepin en allant dans le menu "Runtime -> Run all" ci-haut.

Moreau, D. (2015). When seeing is learning: dynamic and interactive visualizations to teach statistical concepts. Front. Psychol. 6:342. doi: 10.3389/fpsyg.2015.00342


```{code-cell} ipython3
%pip install ipywidgets
%pip install seaborn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import numpy as np
from ipywidgets import interact_manual
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import seaborn as sns
ranfirst = 1

def normal_dist(mu=170, sigma=20):
  plt.figure(figsize=(10, 6))
  x = np.linspace(mu - 10*sigma, mu + 10*sigma, 1000)
  plt.axvline(mu, color='r', linestyle='--')
  plt.plot(x, stats.norm.pdf(x, mu, sigma), label='Distribution de la population')
  plt.title('Distribution de la taille dans la population', fontsize=20)
  plt.ylabel('Densité', fontsize=20)
  plt.xlabel('Taille (cm)', fontsize=20)
  plt.xlim(100, 250)
  plt.tick_params(labelsize=15)
  plt.show()

interactive_plot = interact_manual(normal_dist, mu=(100,240),sigma=(1,40))

```
