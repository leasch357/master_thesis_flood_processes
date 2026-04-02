

clip <- function(x, a_min, a_max) {
  if (x < a_min) {
    return(a_min)
  } else if (x > a_max) {
    return(a_max)
  } else {
    return(x)
  }
}


trimf <- function(x, a, b, c){
  stopifnot(c > a, c > b, b > a)
  pmax(pmin((x - a)/ (b - a), (c - x) / (c - b)), 0)
}



snow_routine <- function(temp, prec, param) {
  
  # """This function simulates a simple, conceptual snow accumulation/melting
  #   model based on a degree day approach.
  # 
  #   Usage:
  #       P, STATES, FLUXES = HBV.snow_routine(param, temp, prec)
  # 
  #   Input:
  #    param = model parameters                              - numpy.ndarray(4, )
  #              1. Ts    = threshold temperature [C]
  #              2. CFMAX = degree day factor [mm/C]
  #              3. CFR   = refreezing factor [-]
  #              4. CWH   = water holding capacity of snow [-]
  #     temp = time series of temperature                    - numpy.ndarray(T, )
  #     prec = time series of precipitation                  - numpy.ndarray(T, )
  # 
  #   Output:
  #        P = time series of simulated flow exiting from    - numpy.ndarray(T, )
  #            the snowpack (as a result of melt-refreezing)
  #            [mm/Dt]
  #   STATES = time series of simulated storages (all in mm) - numpy.ndarray(T,2)
  #            1st column: water content of snowpack
  #                        (snow component)
  #            2nd column: water content of snowpack
  #                        (liquid component)
  #   FLUXES = time series of simulated fluxes (all in mm/Dt)- numpy.ndarray(T,2)
  #            1st column: refreezing
  #            2nd column: snowmelt
  # 
  #   This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin
  #   and T. Wagener at Bristol University (2015).
  #   SAFE is provided without any warranty and for non-commercial use only.
  #   For more details, see the Licence file included in the root directory
  #   of this distribution.
  #   For any comment and feedback, or to discuss a Licence agreement for
  #   commercial use, please contact: francesca.pianosi@bristol.ac.uk
  #   For details on how to cite SAFE in your publication, please see:
  #   https://www.safetoolbox.info """
  
  
  ###########################
  # Recover model parameters:
  ###########################
  Ts    <- param[1]  # Threshold temperature [C]
  CFMAX <- param[2]  # Degree day factor [mm/C]
  CFR   <- param[3]  # Refreezing factor [-]
  CWH   <- param[4]  # Water holding capacity of snow [-]
  
  N <- length(prec)  # Number of time steps
  
  ###########################
  # Initialise variables
  ###########################
  MLT <- numeric(N)  # Snowmelt [mm/Dt]
  RFZ <- numeric(N)  # Refreezing [mm/Dt]
  VS  <- numeric(N)  # Snowpack depth: solid
  VL  <- numeric(N)  # Snowpack depth: liquid
  Q   <- numeric(N)  # Outflow from snowpack
  
  ###########################
  # Snowpack routine
  ###########################
  for (t in 2:N) {
    
    # carry over previous states
    VS[t] <- VS[t - 1]
    VL[t] <- VL[t - 1]
    
    if (temp[t] < Ts) {
      # snowfall
      VS[t] <- VS[t] + prec[t]
      
      # refreezing
      RFZ[t] <- CFR * CFMAX * (Ts - temp[t])
      RFZ[t] <- clip(RFZ[t], 0, VL[t])
      
      VS[t] <- VS[t] + RFZ[t]
      VL[t] <- VL[t] - RFZ[t]
      
    } else {
      # snowmelt
      MLT[t] <- CFMAX * (temp[t] - Ts)
      MLT[t] <- clip(MLT[t], 0, VS[t])
      
      VS[t] <- VS[t] - MLT[t]
      VL[t] <- VL[t] + MLT[t] + prec[t]
    }
    
    # outflow from snowpack
    Q[t] <- VL[t] - (CWH * VS[t])
    Q[t] <- clip(Q[t], 0, VL[t])
    VL[t] <- VL[t] - Q[t]
  }
  
  # Prepare outputs
  STATES <- cbind(VS, VL)
  FLUXES <- cbind(RFZ, MLT)
  
  #   STATES <- cbind(v, vl)
  #   FLUXES <- cbind(rfz, m)    
  #   
  #   robj <- list(P = P, STATES = STATES, FLUXES = FLUXES)
  #   
  #   return(robj)
  #
  robj <- list(P = Q, STATES = STATES, FLUXES = FLUXES)
  return(robj)
  #return(list(Q = Q, STATES = STATES, FLUXES = FLUXES))
}




hbv_sim <- function(P, ept, param, Case){

stopifnot(Case %in% 1:2)

  # """This function simulates the HBV rainfall-runoff model (Seibert, 1997).
  # 
  #   Usage:
  #       Q_sim, STATES, FLUXES = HBV.hbv_sim(param, P, ept, Case, ini)
  # 
  #   Input:
  #     param = vector of model parameters                   - numpy.ndarray(9, )
  #             1. BETA   = Exponential parameter in soil
  #                         routine [-]
  #             2. LP     = evapotranspiration limit [-]
  #             3. FC     = field capacity [mm]
  #             4. PERC   = maximum flux from Upper to Lower
  #                         Zone [mm/Dt]
  #             5. K0     = near surface flow coefficient
  #                         (ratio) [1/Dt]
  #             6. K1     = upper Zone outflow coefficient
  #                         (ratio) [1/Dt]
  #             7. K2     = lower Zone outflow coefficient
  #                         (ratio) [1/Dt]
  #             8. UZL    = near surface flow threshold [mm]
  #             9. MAXBAS = flow routing coefficient [Dt]
  #         P = time series of effective precipitation       - numpy.ndarray(T, )
  #             reaching the ground (i.e. precipitation -
  #              snow accumulation + snowmelt )
  #       ept = time series of potential evapotranspiration  - numpy.ndarray(T, )
  #      Case = flag for preferred path in the Upper Zone    - scalar
  #             dynamics
  #           flag=1 -> Preferred path is runoff
  #           flag=2 -> Preferred path is percolation
  # 
  #   Output:
  #     Q_sim = time series of simulated flow (in mm)        - numpy.ndarray(T, )
  #    STATES = time series of simulated storages            - numpy.ndarray(T,3)
  #             (all in mm)
  #             1: water content of soil (soil moisture)
  #             2. water content of upper reservoir of flow
  #                routing routine
  #             3. water content of lower reservoir of flow  - numpy.ndarray(T,5)
  #                routing routine
  #    FLUXES = time series of simulated fluxes
  #             (all in mm/Dt)
  #             1: actual evapotranspiration
  #             2: recharge (water flux from soil moisture
  #                accounting module to flow routing module)
  #             3: percolation (water flux from upper to
  #                lower reservoir of the  flow routing module)
  #             4: runoff from upper reservoir
  #             5: runoff from lower reservoir
  # 
  #    References:
  # 
  #    Seibert, J.(1997), Estimation of Parameter Uncertainty in the HBV Model,
  #    Nordic Hydrology, 28(4/5), 247-262.
  # 
  #    Comments:
  #    * The Capillary flux (from upper tank to soil moisture accounting module)
  #    is not considered
  #    * The recharge from the soil to the upper zone is considered to be a
  #    faster process than evapotranspiration.
  #    * The preferential path from the upper zone can be modified
  #              - Case 1: interflow is dominant
  #              - Case 2: percolation is dominant
  #    This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin
  #   and T. Wagener at Bristol University (2015).
  #   SAFE is provided without any warranty and for non-commercial use only.
  #   For more details, see the Licence file included in the root directory
  #   of this distribution.
  #   For any comment and feedback, or to discuss a Licence agreement for
  #   commercial use, please contact: francesca.pianosi@bristol.ac.uk
  #   For details on how to cite SAFE in your publication, please see:
  #   https://www.safetoolbox.info"""

# --------------------------
# Recover model parameters:
# --------------------------

BETA <- param[1] # Exponential parameter in soil routine [-]
LP <- param[2] # evapotranspiration limit [-]
FC <- max(.Machine$double.eps, param[3]) # field capacity [mm] cannot be zero

PERC <- param[4] # maximum flux from Upper to Lower Zone [mm/Dt]
K0 <- param[5] # Near surface flow coefficient (ratio) [1/Dt]  
K1 <- param[6] # Upper Zone outflow coefficient (ratio) [1/Dt]  
K2 <- param[7] # Lower Zone outflow coefficient (ratio) [1/Dt]  
UZL <- param[8] # Near surface flow threshold [mm]

MAXBAS <- max(1, floor(param[9])) # Flow routing coefficient [Dt]


N <- length(ept) # number of time samples


# ---------------------
# SOIL MOISTURE ROUTINE
# ---------------------

EA <- numeric(N) # Actual Evapotranspiration [mm/Dt]
#SM <- numeric(N + 1) #  Soil Moisture [mm]
SM <- numeric(N) #  Soil Moisture [mm]
R  <- numeric(N) # Recharge (water flow from Soil to Upper Zone) [mm/Dt]
#UZ <- numeric(N + 1) #  Upper Zone moisture [mm]
#LZ <- numeric(N + 1) #  Lower Zone moisture [mm]
UZ <- numeric(N) #  Upper Zone moisture [mm]
LZ <- numeric(N) #  Lower Zone moisture [mm]
RL <- numeric(N) #  Recharge to the lower zone [mm]
QQ = numeric(N)    # storm flow (quickflow) [mm/Dt]
Q0 <- numeric(N) # Outflow from Upper Zone [mm/Dt]
Q1 <- numeric(N) #  Outflow from Lower Zone [mm/Dt]


#initial stages ini?


for(t in 2:N){
  # start with storages from the previous day
  
  SM[t] = SM[t-1]
  
  UZ[t] = UZ[t-1]
  
  LZ[t] = LZ[t-1]
  
# --------------------------
#    Soil Moisture Dynamics:
# --------------------------

  # precipitation presses onto soil moisture to produce recharge,
  
  # therefore this step comes first  
  
    R[t] <- P[t] * ((SM[t] / FC)^BETA) # Compute the value of the recharge to the 
    R[t] = clip(R[t], 0, P[t])
    
    # we first remove evapotranspiration to reduce excess flow
    EA[t] <- ept[t] * clip(SM[t] / (FC * LP),0, 1) # Compute the evaporation
    EA[t] = clip(EA[t], 0, SM[t])
    
    #update soil moisture storage
    SM[t] = SM[t]+P[t]-R[t]-EA[t]
    
    
    # saturation excess flow represents the overflow of the bucket
    
    QQ[t]  = SM[t] - FC
    QQ[t]  = clip(QQ[t], 0, SM[t])
    SM[t]  = SM[t] - QQ[t]
    
    # make sure that soil moisture stays between zero and field capcity
    # only necessary for numerical stability
    
    SM[t] = clip(SM[t], 0, FC)
    
    
    ####
    #upper zone (we assumed that this process is faster than evaporation)
    #SM_dummy <- max(min(SM[t] + P[t] - R[t], FC), 0) # Compute the water balance 
    # with the value of the recharge  
    #R[t] <- R[t] + max(SM[t] + P[t] - R[t] - FC, 0) + min(SM[t] + P[t] - R[t], 0) #adjust R 
    # by an amount equal to the possible negative SM amount or to the 
    # possible SM amount above FC
    ####
    
    
    #EA[t] <- ept[t] * min(SM_dummy / (FC * LP), 1) # Compute the evaporation
    #SM[t] <- max(min(SM_dummy - EA[t], FC), 0) # Compute the water balance 
    
    #EA[t] <- EA[t] + max(SM_dummy - EA[t] - FC, 0) + min(SM_dummy - EA[t], 0) # adjust EA
    # by an amount equal to the possible negative SM amount or to the 
    # possible SM amount above FC
    
    
    
    
    
    
     
# --------------------
# Upper Zone dynamics: (input first bucket)
# --------------------

    # for the upper zone, recharge is added first, as otherwise there is a
    
    # one day lag between a precipitation event and the runoff response
    
    UZ[t] = UZ[t] + R[t]
    
    if (Case == 1){
    # Case 1: Preferred path = runoff from the upper zone 
  	    #Q0[t] <- max(min(K1 * UZ[t] + K0 * max(UZ[t] - UZL, 0), UZ[t]), 0)     
        #RL[t] <- max(min(UZ[t] - Q0[t], PERC), 0)
      
        # calculate runoff first...
        Q0[t] <- K1 * UZ[t] + K0 * clip(UZ[t] - UZL, 0, UZ[t])
        Q0[t] <- clip(Q0[t], 0, UZ[t])
        UZ[t] <- UZ[t] - Q0[t]
        
        # ...then calculate percolation
        RL[t] <- min(UZ[t], PERC)
        UZ[t] <- UZ[t] - RL[t]
        
      
	} else { #if Case==2
    # Case 2: Preferred path = percolation
        #RL[t] <- max(min(PERC, UZ[t]), 0)
        #Q0[t] <- max(min(K1 * UZ[t] + K0 * max(UZ[t] - UZL, 0), UZ[t] - RL[t]),0)
    	  # calculate percolation first...
    	  RL[t] <- min(UZ[t], PERC)
    	  UZ[t] <- UZ[t] - RL[t]
    	  
    	  # ...then calculate runoff
    	  Q0[t] <- K1 * UZ[t] + K0 * clip(UZ[t] - UZL, 0, UZ[t])
    	  Q0[t] <- clip(Q0[t], 0, UZ[t])
    	  UZ[t] <- UZ[t] - Q0[t]
	  
    }

  #UZ[t+1] <- UZ[t] + R[t] - Q0[t] - RL[t]

# --------------------
# Lower Zone dynamics: 
# --------------------

    #Q1[t] <- max(min(K2 * LZ[t], LZ[t]), 0)
    #LZ[t + 1] <- LZ[t] + RL[t] - Q1[t]
    
    # similarly, recharge to the lower zone is added first, as otherwise
    # there would be an artifical lag again
    
    LZ[t] = LZ[t] + RL[t]
    
    
    # runoff from the lower zone
    
    Q1[t]  = K2*LZ[t]
    
    Q1[t]  = clip(Q1[t], 0, LZ[t])
    
    LZ[t] = LZ[t] - Q1[t]
    
   
}
 
# total runoff is the sum of outflows from upper and lower zone
Q <- QQ + Q0 + Q1 # total outflow (mm/Dt)

# --------------------
# FLOW ROUTING ROUTINE
# --------------------

c <- trimf(1:MAXBAS, 0, (MAXBAS + 1)/2, MAXBAS + 1) #(Seibert,1997)
c <- c / sum(c) # vector of normalized coefficients - (1,MAXBAS)


Q_sim = convolve(Q, c, type = "open")[1:N]

STATES <- cbind(SM, UZ, LZ)
FLUXES <- cbind(EA, R, RL, Q0, Q1)

robj <- list(Q_sim = Q_sim, STATES = STATES, FLUXES = FLUXES)

return(robj)
}

