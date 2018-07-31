# VBM
			CODE DENGUE

(C) Copyright: 2017 V T Romeo Aznar and Hernán G. Solari unless otherwise stated.

Dengue is released under GNU GPL license (see license.txt)

Dengue is a code that simulates the spatial dynamics of an epidemic.

compile as:
gcc -o dengue dengue.c -lm


The main program is dengue.c. (There, there is a short explanation about the system to be simulated)

dengue.inp is an imput file that selects the output. And has the form:
base_file_name_to_use
Number_of_populations Number_of_events
List of populations
List of events
Transitory Simulation "Breeding_sites" (Carring capacity)
where the Transitory and Simulation times are in days (integers).

Some parameters are given at compilation time and are defined in dengue.h

#define LIMfilas 10 /* number of lines in the spatial grid */
#define LIMcolumnas 10  /* Number of columns */
#define PASO 12 /*the day is divided in this number of time steps */
#define REPITE 2 /* Number of repetitions of the same simulation with different
	random series, used to produce statistis */ /* This is a stochastic model. With one repetition, you have one scene. With n repetitions, you have n equivalent scenes.*/
#define TsaveCI /* Transition time (time before the simulation) */ (in this version it is recomended to not use Transition time, just put TsaveCI=0)


In CONSTANTES, are given values for many constants. The human spatial distribution  and carrying capacities are initialized here.
In rates, i) there are rates that change in time (because depend on the population)
          ii) function of the carrying capacity V(N) (Number Of Vectors)
          iii)  Human spatial distribution (taken all the cases the center of the grid as reference, C.centro, see CONSTANTES)
In f-auxiliares file, in the fuction pobla_iniciales, the spatial position of the first case is given. By default it is set in the center of the grid, coincide con el centro de la distribucion espacial.
with C.Inf0 (see CONSTANTES) you can change de amount of first cases. At te moment only one unit has the first infected human, but it is easy to change this in the fuction pobla_iniciales

If you do not have an initial condition file (CIniciales.dat), this it will be created (transition time).

AGREGAR 8-7-2018
* que los rates tienen que estar en 1/dias
* explicar dengue.inp que es cada una de las cosas
* explicar la temperatura
* explicar como ir al caso no dimensional (chequear todas estas cosas)
* Si tenes ganas, pone un listado de como hiciste cada una de las corridas
 a) No espacial, con temperatura constante
 b) No espacial, con temperatura variable
 c) Espacial, con temperatura constante
 d) Espacial, con temperatura variable

* explicar el output: q imprime... 1) la repeticion, 2) la temperatura y dia, 3)


V T Romeo Aznar, Department of Ecology and Evolution, University of Chicago, vromeoaznar at uchicago.edu
-----------------------------------------------------------


