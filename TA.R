
################# 
# UPDATED FALL 2025

Base_ms = 20400.00     # annual (12 month)
Base_PhD = 26400.00   # annual

Tuition_res_9cred = 353.50*9 + 372*9    # annual $6,529
Benefits = 2672+2684   # annual $5,356
Tuition_nonres = 9316 + 9800.5  # annual (for nonresidents) $19,116

MS_yrs = 2
PhD_yrs = 4 

TotVal_res_ms = (Base_ms*(5/6)) *MS_yrs +  (Tuition_res_9cred*MS_yrs)*1 + Benefits*MS_yrs 
TotVal_nonres_ms = (Base_ms*(5/6))*MS_yrs +  (Tuition_res_9cred*MS_yrs)*1 + Benefits*MS_yrs + Tuition_nonres*MS_yrs
TotVal_res_phd = Base_PhD*(5/6)*PhD_yrs +  (Tuition_res_9cred*PhD_yrs)*1 + Benefits*PhD_yrs 
TotVal_nonres_phd = Base_PhD*(5/6)*PhD_yrs +  (Tuition_res_9cred*PhD_yrs)*1 + Benefits*PhD_yrs + Tuition_nonres*PhD_yrs

TotVal_res_ms     # $57,771
TotVal_nonres_ms    # $96,004
TotVal_res_phd      # $145,542
TotVal_nonres_phd    # $212,008


Fees = 5+97+3.5*9+95+60+125+5+18*9 # $580.5



