# TAGS
Transients Analysis of Grounding Systems (TAGS)
Basically, the algorith is the implemantation of the so-called Hibrid Electromagnetic Model applied to grounding grids. It is Implemented in Julia considering Gauss-Legandre integration solution. 
There are two main algorithms --> HEM and MHEM (which is a faster simplification of HEM, or modified HEM)

Frequency Domain
Example 1 --> Compare both HEM and MHEM;
Example 2 --> Compare HEM considering (and not considering) electromagnetic parameters of the soil with the frequency (Model proposed by Alipio and Visacro).
Example 3 --> Compare HEM and MHEM for a horizontal electrode
