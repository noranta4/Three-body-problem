# Three-body-problem
The motion of an asteroid in a two-dimensional space with a star and a planet.

University project • 2014 - Laboratorio di Fisica Computazionale - BSc in Physics, II year

The statement of the problem is in `problem.pdf`, the solution of all points is explained in the document `Antonio_Norelli_solution.pdf`. The code used is `three-body-problem.c`

## Running the program (italian)
Il programma chiede di specificare il punto dell’esonero che si vuole svolgere attraverso il canale standard di input. 
In questo modo vengono settate le corrette condizioni iniziali. Per particolari punti si dovrà specificare anche MB. Ѐ possibile anche scegliere condizioni iniziali a piacere.
Viene chiesto inoltre di specificare il dt e il tempo totale di integrazione. 
Il programma stampa su file sempre lo stesso numero di punti (~10000) affinché il file di output non superi 1MB, 
in questo modo viene preservata la memoria e contenuto notevolmente il tempo di elaborazione (la scrittura su file è
uno dei processi più lenti). Più il tempo di integrazione è lungo più è grande l’intervallo di tempo fra
due punti consecutivi sul file.
I file di output sono nominati “ANpto%.txt” dove al posto di % c’è il particolare punto dell’esonero o lo
0 se si sono settate condizioni iniziali a piacere.

## Examples of trajectories

![hexagon](https://raw.githubusercontent.com/noranta4/Three-body-problem/master/img/hexagon.PNG)
![curls](https://raw.githubusercontent.com/noranta4/Three-body-problem/master/img/curls.PNG)
![strange](https://raw.githubusercontent.com/noranta4/Three-body-problem/master/img/strange.PNG)
![horseshoe](https://raw.githubusercontent.com/noranta4/Three-body-problem/master/img/horseshoe.PNG)
