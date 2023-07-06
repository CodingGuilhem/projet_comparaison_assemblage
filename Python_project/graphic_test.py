import matplotlib.pyplot as plt

# Données pour les lignes
x = [0, 1, 2, 3, 4]  # Valeurs sur l'axe des x (décalage de 0)

y1 = [1, 4, 9, 16, 25]  # Valeurs pour la première ligne
y2 = [0, 2, 4, 6, 8]  # Valeurs pour la deuxième ligne
y3 = [5, 4, 3, 2, 1]  # Valeurs pour la troisième ligne

# Tracer les lignes
plt.plot(x, y1, label='Ligne 1')
plt.plot(x, y2, label='Ligne 2')
plt.plot(x, y3, label='Ligne 3')

# Ajouter une légende
plt.legend()

# Afficher le graphique
plt.show()
