##Genome
SEX_CHROMS = ['chrX', 'chrY']

#Scoring
SIZE_THR = 1 #in Mb
K_THR = 0.11
MODEL_THR = 3
#alpha used to determine weidening threshold, using normal approximation 
FB_ALPHA = 0.1
KSCORE_ALPHA = 0.01


##Chromosome
N_SYMBOL = 'N'
E_SYMBOL = 'E'
U_SYMBOL = 'U'
DEFAULT_N_THRESHOLD = 10
DEFAULT_E_THRESHOLD = 3
N_STR_LEN_THR = 100
HE_Z_THR = 13.8

##Run
WINDOWS_THRESHOLD = 9
SNPS_IN_WINDOW = 1000
WINDOWS_TO_TEST_THRESHOLD = 20
UNIFORMITY_THRESHOLD = 1e-5
LENGTH_THRESHOLD = 10

AI_SENSITIVE_Z = 2
AI_FULL_Z = 2.5
M_Z = 2.5
SINGLE_P_FULL = 0.01
SINGLE_P_SENSITIVE = 0.05

##Segment
MAX_AI_THRESHOLD_FOR_SENSITIVE = 0.1
K_MAX = 1.1

##Testing
COV_INITIAL_SHAPE = 0.14
COV_SHAPE_RANGE = (-2,1)

HE_COV_PERC_BOUNDS = (0.05, 99.0)
HE_VAF_BOUNDS = (0.4,0.6)
HE_FCOV_BOUNDS = (0.01, 0.8)
HE_FN_BOUNDS = (0.2,1.0)
HE_A_BOUNDS = (0.1,0.9)
HE_B_BOUNDS = (1,10)
HE_LERR_BOUNDS = (2,10)


VAF_VAF_BOUNDS = (0.45,0.55)
VAF_N_THR = 100

FB_F_MAX = 1.4
FB_EPS = 1e-4

##Distribution
LENGTH_THRESHOLD = 10
