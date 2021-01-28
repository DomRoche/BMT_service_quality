Data and analysis script for:
Roche DG, Jornod, M, Grutter AS & Bshary R (in prep) Client fish traits underlying variation in service quality in a marine cleaning mutualism

Data collected by DGR & MJ. Please refer to the manuscript for data collection methods and statistical analyses. For questions or to notify the authors if any errors are identified in the data, please contact Dr. Dominique Roche (dominique.roche@mail.mcgill.ca).



### BMT_analysis_script.txt ###

Partially annotated R script for the statistical analysis of "sq_data.csv" and to obtain summary statistics from "raw_data_fastStart_and_mucus.csv" and "raw_data_parasites.csv".


### sq_data.csv ###

Processed data for analysis using the R script BMT_analysis_script.txt. Note that two species of surgeonfishes (Acanthuiridae) were excluded from the analyses (Acanthurus nigrofuscus & Zebrasoma scopas).

columnHeading		description
-------------		-----------
genus			Genus
species			Species
TL			Total length in cm
initiation		Identity of the partner that initiated the interaction : 1=cleaner, 2=client
jolts			Number of jolts during the interaction
end			Identity of the partner that ended the interaction : 1=cleaner, 2=client
endType			Ending of the interaction: 1 = one of the two partners swims away; 2 = client chases the cleaner; 3 = client escapes the cleaner by swimming away at high speed (fast-start);  4 = cleaner switches partner
duration		Duration of the interaction (s)
Desc			Escape distance in a fixed amount of time (average of stage 1 + 2 druration across all fishes) (cm) (Lateral stimulation)
Umax			Maximum velocity within a fixed amount of time (average of stage 1 + 2 druration across all fishes); after smoothing (cm/s) (Lateral stimulation)
Amax			Maximum acceleration within a fixed amount of time (average of stage 1 + 2 druration across all fishes); after smoothing  (cm/s^2) (Lateral stimulation)
turnRate		Stage 1 turning rate as stage 1 turnAngle divided by s1Duration (°/s) (Lateral stimulation)
clientType		Partner choice options: 1=low, 2=intermediate, 3=high (see table A5 in the paper)
mucusMass		Dry mucus mass per unit body surface area (mg/cm^2)
mucusProt		Mucus protein content from a Bradford assay (ug/mg dry weight)
mucusCal		Mucus caloric content (calories/g dry weight) calculated as = -227 + 152 (% carbon) following Platt et al. (1969)
parCm2			Number of non-gnathiid ectoparasites per cm^2 of body surface area
gnaCm2			Number of gnathiid isopods per cm^2 of body surface area



### raw_data_fastStart_and_mucus.csv ###

Raw data for swimming performance (fast-start escape responses) and mucus characteristics for the 13 study species (2 species are removed in the script for Tables A1, A2, and A3). These data were processed for inclusion in sq_data.csv (unfortunately, without using script).

columnHeading		description
-------------		-----------
ID			Fish ID
name			Fish name
genus			Genus
species			Species
mass			Mass (g)
TL			Total length. The straight-line distance from the anterior tip of the snout to the posterior tip of the caudal fin (cm)
SL			Standard length. The straight-line distance from the anterior tip of the snout to the end of the backbone (cm)
BD			Body depth. The maximum straight-line depth of the body (cm)
BW			Body width. The maximum width of the body (cm)
CPW			Caudal peduncle width. The narrowest width measured at the caudal peduncule (cm)
muscleWhite		White muscle mass dissected from one half of the body (g)
muscleRed		White muscle mass dissected from one half of the body (g)
CoM			Center of mass. The straight-line distance from the anterior tip of the snout to the center of mass (cm)
CoM_TL			CoM divided by TL
caudalArea		Caudal fin surface area (cm^2) - for one side of the fin
analArea		Anal fin surface rea  (cm^2) - for one side of the fin
dorsalArea		Dorsal fin surface area (cm^2) - for one side of the fin
dorsalArea2		Second dorsal fin surface area (cm^2) if present - for one side of the fin
pelvicArea		Pelvic fin surface area (cm^2) - for one side of the fin
bodyArea		Body surface area (one side) excluding fins, measured on a flat picture of the fish (cm^2)
pectArea		Pectoral fin surface area (cm^2)  - for one side of the fin
aluArea	Surface 	area of an aluminium sheet laid over the contour of the fish's body (one side, excluding fins) (cm^2)
barbArea		Surface area of one barbel if present (cm^2)
TotCaudalArea		Total caudal fin surface area calculated as 2x caudalArea (cm^2)
TotAnalArea		Total anal fin surface area calculated as 2x analArea (cm^2)
TotDorsalArea		Total dorsal fin surface area calculated as 2x dorsalArea (cm^2
TotDorsalArea2		Total second dorsal fin surface area calculated as 2x dorsalArea2 (cm^2)
TotPelvicArea		Total pelvic fin surface area calculated as 4x pelvicArea (cm^2)
TotPectArea		Total Pectoral fin surface area calculated as 4x pectArea (cm^2)
TotAlulArea		Total aluminium surface area calculated as 2x aluArea (cm^2)
TotBodyArea		Total body surface area calculated as 2x bodyArea (cm^2)
AllAluArea		Sum of all surface areas including the fins and aluminium foil (cm^2)
AllBodArea		Sum of all surface areas including the fins and body measured directly on the pohotograph (cm^2) - should always be lower than allAluArea
relMucusM		Relative mucus mass calculated as mucus mass divided by TotAluArea (mg/cm^2)
mucusProt		Mucus protein content in 1 mg of dried mucus - from Bradford test result (ug/mg)
FStype			Type of fast start escape response resulting from lateral stimulation. S = single-bend C-start in which the tail does not recoil completely after the formation of the C; C-start lacking stage 2. D = double-bend C-start showing a clear full return flip of the tail during stage 2; C-start including both stages 1 and 2. D- = double-bend C-start without full return flip of the tail (without 3rd change in turning direction of the head) or single-bend with muscle contraction during the return of the tail to a straight position. (Lateral stimulation)
latency			Response latency (s). Time between first contact of the stimulus with the water surface and the first head movement (1st contact with water surface = frame 0)
s1Duration		Duration of stage 1 (s). Stage1 = strong unilateral contraction of the body musculature which bends the fish into a C shape
s2Duration		Duration of stage 2 (s). In the case of a single-bend = return of tail until straight (without muscle contraction) - starts with end of stage 1 and end with straightness of the fish (until caudal peduncle, without taking into account caudal fin) - this is not really a stage 2 because there is no stage 2 in single-bends. In the case of a bouble-bend: starts with end of stage 1 and ends with next reversal of the turning direction of the head or the head aligning straight with the body. (Lateral stimulation)
distance		Distance between the fish's CoM and the point of contact of the stimuls with the water as the center of PVC tube in which the stimulus fell (cm) (Lateral stimulation)
stimAngle		Angle of the fish's body relative to the stimulus at the onset of the escape response (0°=front to stimulus, 90°=perpendicular to stimulus, 180°= back to stimulus) (°) (Lateral stimulation)
direction		Direction of the escape relative to the stimulus: -1 = toward stimulus, +1 = away from stimulus (Lateral stimulation)
turnAngle		Stage 1 turning angle as the angle between the straight line joining the tip of the head and the CoM at the onset and end of stage 1 (°) (Lateral stimulation)
turnRate		Stage 1 turning rate as stage 1 turnAngle divided by s1Duration (°/s) (Lateral stimulation)
trajectory		Escape trajectory representing the swimming direction taken by the fish in relation to the stimulus; measured as the angle between the straight line CoM-stimulus at beginning of the response and the line snout-CoM at the end of the response (°) (Lateral stimulation)
turnRadius		Stage 1 turning radius. Area of the circle fitting the best the path made by the center of mass during stage 1; calculated as = (area/pi)^1/2 (cm) (Lateral stimulation)
Desc			Escape distance in a fixed amount of time (average of stage 1 + 2 druration across all fishes) (cm) (Lateral stimulation)
relDesc			Relative escape distance as Desc/TL (Lateral stimulation)
Umax			Maximum velocity within a fixed amount of time (average of stage 1 + 2 druration across all fishes); after smoothing (cm/s) (Lateral stimulation)
Amax			Maximum acceleration within a fixed amount of time (average of stage 1 + 2 druration across all fishes); after smoothing  (cm/s^2) (Lateral stimulation)
F_latency		Frontal stimulation: Response latency as the time between first contact of stimulus with water surface and first head movement (1st contact with water surface = 0). (s)
F_s1Duration		Frontal stimulation: Stage 1 duration (s)
F_s2Duration		Frontal stimulation: Stage 2 duration (s)
F_distance		Frontal stimulation: Distance between the fish's CoM and the point of contact of the stimuls with the water as the center of PVC tube in which the stimulus fell (cm)
F_stimAngle		Frontal stimulation: Angle of the fish's body relative to the stimulus at the onset of the escape response (0°=front to stimulus, 90°=perpendicular to stimulus, 180°= back to stimulus) (°)
F_turnAngle		Frontal stimulation: Direction of the escape relative to the stimulus: -1 = toward stimulus, +1 = away from stimulus
F_turnRate		Frontal stimulation: Stage 1 turning angle as the angle between the straight line joining the tip of the head and the CoM at the onset and end of stage 1 (°)
F_trajectory		Frontal stimulation: Stage 1 turning rate as stage 1 turnAngle divided by s1Duration (°/s)
F_turnRadius		Frontal stimulation: Escape trajectory representing the swimming direction taken by the fish in relation to the stimulus; measured as the angle between the straight line CoM-stimulus at beginning of the response and the line snout-CoM at the end of the response (°)
CHNmass			Mass of the mucus sample used for the CHN analysis.
nitrogen		Percent nitrogen
carbon			Percent carbon
hydrogen		Percent hydrogen
calories		Caloric content of mucus calculaterd as -227+152*carbon (from Platt et al 1969) (Cal/g)



### raw_data_parasites.csv ###

Raw data collected on a sample (1-5) of individuals across the 13 client reef fish species to estimate parasite load (2 species are removed in the script for Table A4). These data were processed for inclusion in sq_data.csv (unfortunately, without using script).

columnHeading		description
-------------		-----------
order			Host fish number in the order in which they were collected
genus			Genus of host
species			Species of host
date			Date of collection day_month_year
site			Site of collection
TL			Total length measured on photograph using ImageJ (cm)
SL			Standard length measured on photograph using ImageJ (cm)
bodySurface		Body surface area (one side of the body) excluding fins measured on photograph using ImageJ (cm^2)
regLSlope		Slope of the linear regression between BodySurface and total surface area (frommucus data)
regLIntercept		Intercept of the linear regression between BodySurface and total surface area (from mucus data)
regLRsquared		R squared of the linear regression between BodySurface and total surface area (from mucus data)
totalSurface		Total surface area calculated using the linear regression between BodySurface and total surface area (mucus data)
ano1			Parasite count: Anoplodiscus-like monogenean species 1 (smaller than Ano2)
ano2			Parasite count: Anoplodiscus-like monogenean species 1 (larger than Ano1)
bomo			Parasite count: Bomolochidae sp.
cali			Parasite count: Caligidae sp.
copLarv			Parasite count: copepod larva (smaller than Bomo and Cali and mainly present in H. melapterus)
dactly			Parasite count: Dactylogiridae sp.
flatMono		Parasite count: flat monogenean with round head and prominent hooks
gna			Parasite count: gnathiid isopod
nema			Parasite count: nematode (probably from feces in collection bag) should not be counted as an ectoparasite
roundMono		Parasite count: monogenean with round (pinkish) body and round head with prominent hooks (similar to flat mono but body is much rounder and different colour, more opaque)
unIdCop			Parasite count: unidentified parasitic copepod
unIdMono		Parasite count: unidentified monogenean
unID			Parasite count: unidentified parasite
notes			Notes
totParasite		Total number of ectoparasites
totParasiteCm2		Number of ectoparasites per cm^2 of body surface area
otherParasiteCm2	Number of non-gnathiid ectoparasites per cm^2 of body surface area
gnaCm2			Number of gnathiid isopods per cm^2 of body surface area



### raw_data_field_observations.csv ###

Raw data of cleaning interactions observed in the field for the 13 client reef fish species (2 extra species are included). These data were processed for inclusion in sq_data.csv (unfortunately, without using script).

columnHeading		description
-------------		-----------
genus			Genus
species			Species
TL			Total length as the straight-line distance from the anterior tip of the snout to the posterior tip of the caudal fin (cm)
initiation		Identity of the partner that initiated the interaction : 1=cleaner, 2=client
jolts			Number of jolts during the interaction
end			Identity of the partner that ended the interaction : 1=cleaner, 2=client
endType			Ending of the interaction: 1 = one of the two partners swims away; 2 = client chases the cleaner; 3 = client escapes the cleaner by swimming away at high speed (fast-start);  4 = cleaner switches partner
duration		Duration of the interaction (s)
date			Date of observation
time			Approximate time of observation
site			Site of observation
method			Method of observation : snorkling or diving
observer		DR = Dominique Roche, MJ = Maïwenn Jornod
jolts100sec		Number of jolts (jolts) divided by the duration of the interaction (duration) multiplied by 100


### References

Platt, T., Brawn, V. M., & Irwin, B. (1969). Caloric and carbon equivalents of zooplankton biomass. Journal of the Fisheries Board of Canada, 26(9), 2345-2349. 
