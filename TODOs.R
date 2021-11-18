#For function storage/ developement
#TODOs

#1) Prevent mother-son mating: probably some mother-son mating in nature as well - so that's fine

# Functions
#1) Write function reset events
#2) Create a function to sample the location for the swarm: later
#3) Create a function to sample locations within a radius: later
#4) Write a setPheno function for the colony: later
#5) Revise last event slot (whether we use/need it)
#7) THink about how to handle the colour
#8) Call setPheno after swarm, split and supersede, crossColony, createColony (follow the AlphaSimR logic of SP parameter vs setPheno)
#9) Assign csd locus - look Laura's code
#10) Make csd function in terms of inbreeding
#11) Add names to colonies in the Colonies (name = colony id)
#12) Consider adding a combine colony function (puts workers from weak into the strong colony)
#16) Think about use = "rand/something" in the level 3 functions (no need for ID then)
#14) All functions should test the class and throw a stop (colony, colonies)
#17) Only set the id when we have the queen!
#18) All add functions should call create functions
#19) Do we just remove provided qeens and fathers in colony (level2) functions? 
#20) CHange value 10 in createMultipleMatedColonies (put in the arguments?)
#21) Revise collapse slot and functions (should we just remove them?)
#22) Check the sequence of events (before turning production on) - think about event table (3x3x3 matrix for 3 events)

# Text
#1) Think about providing informative messages for the functions: Laura
#2) Think of a good names for the swarmed colony (the one that stay)

# Think
#1) Think about removing workers and drones in "instantaneous" functions (opposite to adding them)
#2) Think about replacing phenotypes in the swarm/supersede/split


#TODO for script
#S1) Distribute supersedure/swarming events throughout the year/seasons
#S2) Distribute colony losses events throughout the year/seasons

