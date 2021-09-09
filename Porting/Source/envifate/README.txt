EnviRisk Qgis plugin

Tutto il pacchetto va copiato/clonato in 
/home/user/.qgis2/python/plugin per gli utenti Linux
/User/nomeutente/.qgis2/python/plugin per gli utenti Mac
C:\Users\(user)\.qgis2\python\plugin per gli utenti Windows

Il plugin è strutturato nella classica struttura a cartelle standard dei plugin Qgis, con in più una cartella
 library, che contiene tutti i moduli di analisi di rischio che rappresentano il cuore dell'applicazione.

In particolare ad oggi nella cartella library sono presenti i seguenti files:
daf.py -> calcola l'attenuazione laterale in falda
leaching.py -> calcola l'infiltrazione in falda delle sostanze inquinanti
river.py -> simula la dinamica 1d di concentrazione di un inquinante in un corso d'acqua
lake.py -> modello 2d per la simulazione della dinamica di un inquinante in un lago
plume.py -> calcola la concentrazione e il livello del pennacchio inquinante da una sorgente tipo "camino"
functions.py -> raccoglie tutte le funzioni generali che vengono riutilizzate dai diversi moduli
substance.db -> database in sqlite che contiene tutti i parametri che vengono utilizzati dai moduli

Ogni modulo corrispondente al file.py viene "caricato" dal plugin a seconda del processo di analisi 
selezionato (ad oggi dentro il plugin funziona esclusivamente il modulo leaching). 
La natura "a libreria" permette tuttavia di lanciare da terminale oppure da script 
o applicazione esterna ogni singolo modulo.

Da terminale è possibile conoscere i parametri necessari per ciascun plugin eseguendo lo stesso in modalità
 "aiuto" ovvero posponendo il parametro "-h". Ad esempio
python leaching.py -h 
mostra la descrizione di tutti i parametri necessari per effettuare il calcolo.

Autori: Francesco Geri, Marco Ciolli
Consulenti: Paolo Zatelli, Oscar Cainelli

Licenza: GPL v.3