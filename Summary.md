# Summary

## Toy Image - Cellule
Immagine di prova di alcuni _acini_ senza nessuna relazione spaziale tra di loro e zoom di una regione:
  - Ogni poligono rappresenta una _cellula_ con il proprio _nucleo_.
  - Le regioni più scure rappresentano le pareti.
  - Le regioni grigie rappresentano le regioni interne.
  - La regione bianca rappresenta lo spazio tra i diversi acini.

  <img src="/Segmentation/Zoom/full_zoom.png" width="450"/> <img src="/Segmentation/Zoom/bw_nuclei_bound.png" width="400"/>

#### Segmentazione
Basandosi, per ora, sulle diverse tonalità di colore si riescono a distinguere le diverse regioni, a cui si può assegnare una sorta di identià:

<img src="/Segmentation/Zoom/bounds.png" width="400"/> <img src="/Segmentation/Zoom/lumes.png" width="400"/>

<img src="/Segmentation/Zoom/cells.png" width="400"/> <img src="/Segmentation/Zoom/nuclei.png" width="400"/>

## Distribuzione spaziale
Per cercare di ricreare l'organizzazione spaziale del tessuto si è pensato di usare il seguente pattern:

<img src="/Ramification/ramification_6.png" />

Che può essere ramificato a piacere, ottenendo per esempio un disegno del genere:

<img src="/Ramification/ramification_13_wo_noise.png" width=""/>

Inserendo una leggera randomizzazione sulla rimificazione del disegno:

<img src="/Ramification/ramification_13.png" width=""/>
