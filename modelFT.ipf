function modelFT(SCS, SCD, Name, PARAM)
//Function for two-compartment model:  Inputs are SCS and SCD, outputs waves that start with the prefix
// supplied by the string "Name" and include the somatic voltage (Name_ES), dendritic voltage(Name_ED), 
//potassium conductance (Name_GKS), and delta functions whenever spikes occur (Name_SP)
// model parameter values are supplied by the wave PARAM (sample values are: (0.1,5,-15,0.036,55,0.06,0.024,0.24,0.018,0.06,12,28000)
wave SCS // Input wave.  Current injected into soma simulates a pulse of current
wave SCD // Input wave.  Current injected into dendrite
string Name // Input string.  Defines prefix of output waves
wave PARAM // Input wave.  Wave contains all other input parameters
variable step = PARAM[0] // Input variable.  Integration step size in ms.
variable fine = PARAM[1] // Input variable.  Number of fine steps per step above.
variable ek = PARAM[2] // Input variable.  Potassium equilibrium potential relative to rest in mV.
variable b = PARAM[3] // Input variable.  Potassium conductance amplitude due to AP in uS
variable tgk = PARAM[4] // Input variable.  Potassium conductance decay time constant in ms.
variable cms = PARAM[5] // Input variable.  Soma capacitance in nF.
variable gls = PARAM[6] // Input variable.  Soma conductance in uS.
variable cmd = PARAM[7] //Capacitance of dendritic compartment
variable gld = PARAM[8] // Input variable.  Dendrite conductance in uS.
variable gc = PARAM[9] // Input variable.  Coupling conductance in uS.
variable tho = PARAM[10] // Input variable.  Minimum threshold of soma in mV relative to rest.
variable ltstop = PARAM[11] // Input variable.  Simulation time in ms.
variable refract = 2.0 //absolute refractory period in ms

variable i,j // Index variables for nested loops.
variable length = ltstop / step //Size of output waves.
variable esn, edn // Temporary variables to hold es & ed during fine integration.
variable p, q, r, s // Temporary variables to for integrations
variable dtfine = step/fine // Time (in ms) of fine steps.
variable refractory_counter=0 // Counter used to keep track of refractory periods.

string Es_wavename; Es_wavename = Name + "_ES"; make /o/d/n=(length) $Es_wavename
string Sp_wavename; Sp_wavename = Name + "_SP"; make /o/d/n=(length) $Sp_wavename
string Gk_wavename; Gk_wavename = Name + "_GKS";make /o/d/n=(length) $Gk_wavename
string Ed_wavename; Ed_wavename = Name + "_ED"; make /o/d/n=(length) $Ed_wavename

wave ES = $Es_wavename // Output wave.  Membrane potential of soma
wave SP = $Sp_wavename // Output wave.  Spiking variable of soma
wave GKS = $Gk_wavename // Output wave.  Potassium conductance of soma
wave ED = $Ed_wavename // Output wave.  Membrane potential of dendrites

ES[0] = 0
SP[0] = 0
GKS[0] = 0
ED[0] = 0

i = 1
do
j = 0
esn = ES[i-1]
edn = ED[i-1]

if((esn > tho) %& (refractory_counter <= 0))
SP[i] = 1;
refractory_counter = refract/step
else
SP[i] = 0
endif

if(refractory_counter <= 0) // that is if a spike has occurred more than 0.5 ms ago

//Updating the differential equations is based on the exponential integration scheme described in MacGregor, R.J. Neural 
// and Brain Modeling. Academic Press, NY, 1987.  The basic idea is as follows:
// if a differential equation can be put in the following form: dV/dt = -A*V +B then the value of V can be calculated for each
// successive time step (of length dt) as follows: V[i+1]=V[i]*exp(-A*dt) +(B/A)*(1-exp(-A*dt).  For example, the
// differential equation for somatic voltage is as follows:  desn/dt = (1/cms)*(-gls*esn -gc*(esn - edn) - GKS*(esn-ek) +SCS)
// rearranging the equation we get:  desn/dt = -((gls + gc +GKS)/cms)*esn + ((SCS +gc*edn +GKS*ek)/cms)
// if we let p = (gls+gc+GKS[i-1])/cms and q = (SCS[i-1] + gc*edn+GKS[i-1]*ek)/cms, then the equation is updated as described 
// below.  A similar equation is obtained for updating the dendritic voltage.  Because current flow between the two compartments is
// relatively fast, we update esn and edn several times for each major time step.
do
p = (gls+gc+GKS[i-1])/cms
q = (SCS[i-1] + gc*edn+GKS[i-1]*ek)/cms

r = (gld+gc)/cmd
s = (SCD[i-1]+gc*esn)/cmd

esn = esn*exp(-p*dtfine) + (q/p)*(1-exp(-p*dtfine))
edn = edn*exp(-r*dtfine) + (s/r)*(1-exp(-r*dtfine))

j += 1
while(j < fine)
else
refractory_counter -= 1
//  if a spike has occurred less than 0.5 ms ago, we temporarily ignore the stimulus current and the potassium current
// this may be different in the eventual NEURON version of the model, put should lead to relatively small differences in behavior
do
p = (gls+gc)/cms
q = ( gc*edn)/cms

r = (gld+gc)/cmd
s = (gc*esn)/cmd

esn = esn*exp(-p*dtfine) + (q/p)*(1-exp(-p*dtfine))
edn = edn*exp(-r*dtfine) + (s/r)*(1-exp(-r*dtfine))

j += 1
while(j < fine)
endif
ES[i] = esn
ED[i] = edn
GKS[i] = GKS[i-1]*exp(-step/tgk) + b*SP[i]
i+= 1
while(i<length)
end
