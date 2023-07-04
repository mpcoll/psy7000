{
 "cells": [
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {
    "id": "esC19go4jvVs"
   },
   "source": [
    "# Cours 1 - Introduction et fondements de l'expérience\n",
    "Créé par Michel-Pierre Coll - 2022, michel-pierre.coll@psy.ulaval.ca\n",
    "\n",
    "\n",
    "## 0. Instructions - calepin interactifs\n",
    "\n",
    "Les simulations, la visualisation et l'interaction sont d'excellentes façon de comprendre certains concepts statistiques (Moreau, 2015). Au cours de la session nous utiliserons régulièrement des simulations et visualisations pour illuster certains concepts. \n",
    "\n",
    "Ces simulations sont effectuées en utilisant des documents comme celui-ci (\"calepins\" ou *notebooks*) qui contiennent des explications et des \"cellules\" qui permettent d'exécuter le langage de programmation *Python*. Dans la majorité des cas, le code sera caché. Il vous sera parfois demandé de changer quelques paramètres dans le code pour évaluer leur effet sur les résultats à l'aide de boutons. **LA COMPRÉHENSION DU CODE DANS LES CELLULES N'EST PAS NÉCESSAIRE POUR LE COURS.** Si vous aimeriez comprendre comment les simulations sont effectuées, vous pouvez cliquer sur \"Show code\".\n",
    "\n",
    "\n",
    "Pour débuter, exécutez le calepin en allant dans le menu \"Runtime -> Run all\" ci-haut.\n",
    "\n",
    "Moreau, D. (2015). When seeing is learning: dynamic and interactive visualizations to teach statistical concepts. Front. Psychol. 6:342. doi: 10.3389/fpsyg.2015.00342\n",
    "\n"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "```{tip} ipython3\n",
    "test\n",
    "\n",
    "```"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "cellView": "form",
    "id": "hxCUg3lC7ylq",
    "tags": [
     "hide-input"
    ]
   },
   "outputs": [],
   "source": [
    ":tags: [hide-input]\n",
    "%pip install ipywidgets\n",
    "%pip install seaborn\n",
    "\n",
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "from scipy import stats\n",
    "import numpy as np\n",
    "import seaborn as sns\n",
    "from ipywidgets import interact_manual\n",
    "%config InlineBackend.figure_format='retina'\n",
    "\n",
    "\n",
    "def normal_dist(mu=170, sigma=20):\n",
    "  plt.figure(figsize=(10, 6))\n",
    "  x = np.linspace(mu - 10*sigma, mu + 10*sigma, 1000)\n",
    "  plt.axvline(mu, color='r', linestyle='--')\n",
    "  plt.plot(x, stats.norm.pdf(x, mu, sigma), label='Distribution de la population')\n",
    "  plt.title('Distribution de la taille dans la population', fontsize=20)\n",
    "  plt.ylabel('Densité', fontsize=20)\n",
    "  plt.xlabel('Taille (cm)', fontsize=20)\n",
    "  plt.xlim(100, 250)\n",
    "  plt.tick_params(labelsize=15)\n",
    "  plt.show()\n",
    "\n",
    "interactive_plot = interact_manual(normal_dist, mu=(100,240),sigma=(1,40))\n",
    "\n"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {
    "id": "G-FT_6BhJqbV"
   },
   "source": [
    "Nous utiliserons souvent des figures interactives. Par exemple, dans la cellule ci-dessous, vous pouvez modifier la moyenne et l'écart-type de la distribution normale de la taille dans la population et voir l'effet de ces changements sur la forme de la distribution.\n",
    "\n",
    "**mu:** Moyenne de la distribution dans la population\n",
    "\n",
    "**sigma:** Écart-type de la distribution dans la population.\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "cellView": "form",
    "id": "4LGMHuXmASjM",
    "tags": [
     "remove-input"
    ]
   },
   "outputs": [],
   "source": [
    "\n",
    "def normal_dist(mu=170, sigma=20):\n",
    "  plt.figure(figsize=(10, 6))\n",
    "  x = np.linspace(mu - 10*sigma, mu + 10*sigma, 1000)\n",
    "  plt.axvline(mu, color='r', linestyle='--')\n",
    "  plt.plot(x, stats.norm.pdf(x, mu, sigma), label='Distribution de la population')\n",
    "  plt.title('Distribution de la taille dans la population', fontsize=20)\n",
    "  plt.ylabel('Densité', fontsize=20)\n",
    "  plt.xlabel('Taille (cm)', fontsize=20)\n",
    "  plt.xlim(100, 250)\n",
    "  plt.tick_params(labelsize=15)\n",
    "  plt.show()\n",
    "\n",
    "interactive_plot = interact_manual(normal_dist, mu=(100,240),sigma=(1,40))\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "#@title _\n",
    "\n",
    "def normal_dist_sample(mu=170, sigma=20, n=20, nouveau=0):\n",
    "    groupe1 = np.random.normal(loc=mu, scale=sigma, size=n)\n",
    "    fig, axes= plt.subplots(1, 2, figsize=(20, 6))\n",
    "\n",
    "    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 1000)\n",
    "    y1 = [stats.norm.pdf(x, mu, sigma)[np.argmin(np.abs(x - v))] for v in groupe1]\n",
    "    axes[0].axvline(mu, color='r', linestyle='--')\n",
    "\n",
    "    axes[0].plot(x, stats.norm.pdf(x, mu, sigma), label='Distribution de la population')\n",
    "    axes[0].set_title('1A- Distribution de la taille dans la population', fontsize=20)\n",
    "    axes[0].scatter(groupe1, y1, label=('Dernier échantillon aléatoire (N=' +  str(n) + ')') , s=50, color='k')\n",
    "    axes[0].set_ylabel('Densité', fontsize=20)\n",
    "    axes[0].set_xlabel('Taille (cm)', fontsize=20)\n",
    "    axes[0].set_xlim(mu-3*sigma, mu+3*sigma)\n",
    "    axes[0].tick_params(labelsize=15)\n",
    "    axes[0].legend(fontsize=12)\n",
    "\n",
    "    sns.histplot(groupe1, kde=False, edgecolor='k', ax=axes[1])\n",
    "    axes[1].axvline(np.mean(groupe1), linestyle='--', color='r', label='Moyenne: (N=' + str(n) + '): ' + str(np.round(np.mean(groupe1), 2)))\n",
    "    axes[1].axvline(np.mean(groupe1), linestyle='--', color='r', alpha=0, label='ÉT: (N=' + str(n) + '): ' + str(np.round(np.std(groupe1), 2)))\n",
    "\n",
    "    # axes[1].axvline(np.mean(all_means)-np.std(groupe1), linestyle='--', color='g', label='Écart-type: ' + str(np.round(np.std(groupe1), 2)))\n",
    "    # axes[1].axvline(np.mean(all_means)+np.std(groupe1), linestyle='--', color='g', label='Erreur standard: ' + str(np.round(np.std(groupe1)/np.sqrt(N_groupe1), 2)))\n",
    "    axes[1].set_title(\"1B - Distribution de l'échantillon aléatoire\", fontsize=20)\n",
    "    axes[1].set_ylabel('Fréquence', fontsize=20)\n",
    "    axes[1].set_xlabel('Taille (cm)', fontsize=20)\n",
    "    axes[1].set_xlim(mu-3*sigma, mu+3*sigma)\n",
    "\n",
    "    axes[1].legend(fontsize=12)\n",
    "    axes[1].tick_params(labelsize=15)\n",
    "    plt.show()\n",
    "\n",
    "\n",
    "interactive_plot = interact_manual(normal_dist_sample, mu=(100,240),sigma=(1,500), n=(1, 5000), nouveau=(0,1))\n",
    "interactive_plot.widget.children[4].description = 'Exécuter!'\n",
    "\n",
    "\n"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## TEst"
   ]
  }
 ],
 "metadata": {
  "colab": {
   "collapsed_sections": [],
   "provenance": []
  },
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.13"
  },
  "vscode": {
   "interpreter": {
    "hash": "5820ac161b8e279ce7635811cf4400d19b4a14f7311591c3c8530b3b01883bed"
   }
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
