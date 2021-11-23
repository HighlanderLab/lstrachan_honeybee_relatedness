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
#23) Write a function set location (setLocation(x) x = Colony|Colonies))
#24) Add "crossPlan" argument to the crossCOlonies, default = "rand"
     # get the fathers from rand_location --> some sort of a plan
#25) Store inheritance, selection and production criteria in the Colony
#26) Make drones properly haploid!
#27) Think about how to name the traits in ALphaSimR and AlphaSimRBee

#Pheno related TODOs
#28) Create getPhenoColony --> returns a vector and getPhenoCOlonies --> returns a matrix
#29) Write function getEventsColony and getEventsColonies
#30) THink about creating getPheno and getEvents (whether they work on colony or colonies)
#31) getGVColony: input is gonna be a colony, output will be a list of queen|workers|drones slot and
# the input is gonna be a FUN --> we can either provide mean | identity 
# maybe the function average is actually applied in another functions
#32) getBvColony, getDdColony
#33) Write a document of the potential to get queen, worker and drone effect for all individuals
#34) Think about indirect (social) genetic effects 
#35) THis about pairing production|events with pheno

# Text
#1) Think about providing informative messages for the functions: Laura
#2) Think of a good names for the swarmed colony (the one that stay)

# Think
#1) Think about removing workers and drones in "instantaneous" functions (opposite to adding them)
#2) Think about replacing phenotypes in the swarm/supersede/split


#TODO for script
#S1) Distribute supersedure/swarming events throughout the year/seasons
#S2) Distribute colony losses events throughout the year/seasons

