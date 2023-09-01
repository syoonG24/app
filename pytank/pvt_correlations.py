import numpy as np
def Bo_bw(P, T, salinity, unit=1):
    """This function calculate the Bo of methane dissolved brine
    The user must select UNITS 1 to Field units and 2 to International System

    IF UNITS 1 is selected introduce:
        | Pressure (P) in PSI
        | Temperature (T)2 in ºF
        | Salinity in ppm
        | Bo returns in bbl/STB

    IF UNITS 2 is selected introduce:
        Pressure (P) in MPa
        Temperature (T) in ºC
        Salinity in ppm
        Bo returns in cm3/cm3"""

    # Step1: Calculate the derivatives of Duan

    if unit == 1:
        P_array = np.array(P)
        P_array = P_array * (6.89475729 * 10 ** -3.0)
        T_array = np.array(T)
        T_array = (T_array - 32) * 5 / 9
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    else:
        P_array = np.array(P)
        T_array = np.array(T)
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    T_array = T_array + 273
    # Coefficients Table 4_11
    UCH4c6 = 7.6985890 * 10 ** (-2)
    UCH4c7 = -5.0253331 * 10 ** (-5)
    UCH4c8 = -30.092013
    UCH4c9 = 4.8468502 * 10 ** (3)

    CH4c6 = 3.92 * 10 ** (-4)
    CH4c10 = -1.97 * 10 ** (-6)

    # Ec 4.19
    dch4dp = UCH4c6 + UCH4c7 * T_array + UCH4c8 / T_array + UCH4c9 / (T_array) ** 2
    # Ec 4.20
    dlamdp = CH4c6 + 2 * CH4c10 * P_array
    # Ec 4.22 Partial molar volume of methane in brine
    VMCH4B = 8.314467 * T_array * (dch4dp + 2 * m * dlamdp + m ** 2 * 0)

    # Step2: Calculate the density of methane-free brine water

    T_array = T_array - 273

    # Coefficients Table 4_6 and Ec 4.1
    a1_T_70 = -0.127213
    a2_T_70 = 0.645486
    a3_T_70 = 1.03265
    a4_T_70 = -0.070291
    a5_T_70 = 0.639589

    a_T_70 = (a1_T_70 * (T_array / 100) ** 2 + a2_T_70 * (T_array / 100) + a3_T_70) / (
                a4_T_70 * (T_array / 100) ** 2 + a5_T_70 * (T_array / 100) + 1)

    a1_EW = 4.221
    a2_EW = -3.478
    a3_EW = 6.221
    a4_EW = 0.5182
    a5_EW = -0.4405

    a_T_EW = (a1_EW * (T_array / 100) ** 2 + a2_EW * (T_array / 100) + a3_EW) / (
                a4_EW * (T_array / 100) ** 2 + a5_EW * (T_array / 100) + 1)

    a1_FW = -11.403
    a2_FW = 29.932
    a3_FW = 27.952
    a4_FW = 0.20684
    a5_FW = 0.3768

    a_T_FW = (a1_FW * (T_array / 100) ** 2 + a2_FW * (T_array / 100) + a3_FW) / (
                a4_FW * (T_array / 100) ** 2 + a5_FW * (T_array / 100) + 1)

    # Coefficients Table 4_7 and Ec 4.1

    a1_DM2 = -1.1149 * 10 ** (-4)
    a2_DM2 = 1.7105 * 10 ** (-4)
    a3_DM2 = -4.3766 * 10 ** (-4)
    a4_DM2 = 0
    a5_DM2 = 0

    a_T_DM2 = (a1_DM2 * (T_array / 100) ** 2 + a2_DM2 * (T_array / 100) + a3_DM2) / (
            a4_DM2 * (T_array / 100) ** 2 + a5_DM2 * (T_array / 100) + 1)

    a1_DM3_2 = -8.878 * 10 ** (-4)
    a2_DM3_2 = -1.388 * 10 ** (-4)
    a3_DM3_2 = -2.96318 * 10 ** (-3)
    a4_DM3_2 = 0
    a5_DM3_2 = 0.51103

    a_T_DM3_2 = (a1_DM3_2 * (T_array / 100) ** 2 + a2_DM3_2 * (T_array / 100) + a3_DM3_2) / (
                a4_DM3_2 * (T_array / 100) ** 2 + a5_DM3_2 * (T_array / 100) + 1)

    a1_DM1 = 2.1466 * 10 ** (-3)
    a2_DM1 = 1.2427 * 10 ** (-2)
    a3_DM1 = 4.2648 * 10 ** (-2)
    a4_DM1 = -8.1009 * 10 ** (-2)
    a5_DM1 = 0.525417

    a_T_DM1 = (a1_DM1 * (T_array / 100) ** 2 + a2_DM1 * (T_array / 100) + a3_DM1) / (
                a4_DM1 * (T_array / 100) ** 2 + a5_DM1 * (T_array / 100) + 1)

    a1_DM1_2 = 2.366 * 10 ** (-4)
    a2_DM1_2 = -3.636 * 10 ** (-4)
    a3_DM1_2 = -2.278 * 10 ** (-4)
    a4_DM1_2 = 0
    a5_DM1_2 = 0

    a_T_DM1_2 = (a1_DM1_2 * (T_array / 100) ** 2 + a2_DM1_2 * (T_array / 100) + a3_DM1_2) / (
                a4_DM1_2 * (T_array / 100) ** 2 + a5_DM1_2 * (T_array / 100) + 1)

    # Ec 4.6 Density of brine water at (70MPa and T)
    dbw_T_70 = a_T_70 + a_T_DM2 * 1 * m ** 2 + a_T_DM3_2 * 1 * m ** (3 / 2) + a_T_DM1 * 1 * m + a_T_DM1_2 * 1 * m ** (
                1 / 2)

    # Coefficients Table 4_8 Ec 4.1

    a1_EM = 0
    a2_EM = 0
    a3_EM = 0.1249
    a4_EM = 0
    a5_EM = 0

    a_T_EM = (a1_EM * (T_array / 100) ** 2 + a2_EM * (T_array / 100) + a3_EM) / (
            a4_EM * (T_array / 100) ** 2 + a5_EM * (T_array / 100) + 1)

    a1_FM3_2 = -0.617
    a2_FM3_2 = -0.747
    a3_FM3_2 = -0.4339
    a4_FM3_2 = 0
    a5_FM3_2 = 10.26

    a_T_FM3_2 = (a1_FM3_2 * (T_array / 100) ** 2 + a2_FM3_2 * (T_array / 100) + a3_FM3_2) / (
            a4_FM3_2 * (T_array / 100) ** 2 + a5_FM3_2 * (T_array / 100) + 1)

    a1_FM1 = 0
    a2_FM1 = 9.917
    a3_FM1 = 5.1128
    a4_FM1 = 0
    a5_FM1 = 3.892

    a_T_FM1 = (a1_FM1 * (T_array / 100) ** 2 + a2_FM1 * (T_array / 100) + a3_FM1) / (
                a4_FM1 * (T_array / 100) ** 2 + a5_FM1 * (T_array / 100) + 1)

    a1_FM1_2 = 0.0365
    a2_FM1_2 = -0.0369
    a3_FM1_2 = 0
    a4_FM1_2 = 0
    a5_FM1_2 = 0

    a_T_FM1_2 = (a1_FM1_2 * (T_array / 100) ** 2 + a2_FM1_2 * (T_array / 100) + a3_FM1_2) / (
                a4_FM1_2 * (T_array / 100) ** 2 + a5_FM1_2 * (T_array / 100) + 1)

    Eb_T_M = a_T_EW + a_T_EM * m

    FB_T_M = a_T_FW + a_T_FM3_2 * m ** (3 / 2) + a_T_FM1 * m + a_T_FM1_2 * m ** (1 / 2)

    Ib_T_70_m = (1 / Eb_T_M) * np.log(np.abs(Eb_T_M + FB_T_M))

    Ib_T_P_m = (1 / Eb_T_M) * np.log(np.abs(Eb_T_M * (P_array / 70) + FB_T_M))

    # Density of methane-free brine wate
    denbwnogas = (dbw_T_70 * np.exp((Ib_T_P_m - Ib_T_70_m)))

    # Step3: Calculate the vapor pressure of pure water

    T_array = T_array + 273

    # Coefficients Table 4_9
    a1Po = -7.85951783
    a2Po = 1.84408259
    a3Po = -11.7866497
    a4Po = 22.6807411
    a5Po = -15.9618719
    a6Po = 1.80122502

    # Ec 4.13
    vad = 1 - T_array / 674.096  ##Adimensional T EC.4.14
    # Ec 4.14
    Po = 22.064 * (np.exp((674.096 / T_array) * (
                a1Po * vad + a2Po * vad ** 1.5 + a3Po * vad ** 3 + a4Po * vad ** 3.5 + a5Po * vad ** 4 + a6Po * vad ** 7.5)))

    # Step4: Calculate the solubility of methane in pure water at (P,T)

    T_array = T_array - 273

    # Coefficients Table 4_10 EC 4.1
    a1A_T = 0
    a2A_T = -0.004462
    a3A_T = -0.06763
    a4A_T = 0
    a5A_T = 0

    A_T = (a1A_T * (T_array / 100) ** 2 + a2A_T * (T_array / 100) + a3A_T) / (
                a4A_T * (T_array / 100) ** 2 + a5A_T * (T_array / 100) + 1)

    a1B_T = -0.03602
    a2B_T = 0.18917
    a3B_T = 0.97242
    a4B_T = 0
    a5B_T = 0

    B_T = (a1B_T * (T_array / 100) ** 2 + a2B_T * (T_array / 100) + a3B_T) / (
                a4B_T * (T_array / 100) ** 2 + a5B_T * (T_array / 100) + 1)

    a1C_T = 0.6855
    a2C_T = -3.1992
    a3C_T = -3.7968
    a4C_T = 0.07711
    a5C_T = 0.2229

    C_T = (a1C_T * (T_array / 100) ** 2 + a2C_T * (T_array / 100) + a3C_T) / (
                a4C_T * (T_array / 100) ** 2 + a5C_T * (T_array / 100) + 1)

    # Ec 4.13 (Solubility of methane in pure water)
    mCH4pw = np.exp(A_T * (np.log(P_array - Po)) ** 2 + B_T * np.log(P_array - Po) + C_T)  ##Ec 4.15

    T_array = T_array + 273

    # Step5: Calculate the solubility of methane in brine water at (P,T)

    # Coefficients Table 4_11
    lambC1 = -0.80898
    lambC2 = 1.0827 * 10 ** (-3)
    lambC3 = 183.85
    lambC6 = 3.924 * 10 ** (-4)
    lambC10 = -1.97 * 10 ** (-6)

    # Ec 4.16
    lambCH4 = lambC1 + lambC2 * T_array + lambC3 / T_array + lambC6 * P_array + lambC10 * P_array ** 2  # EC4.16

    # Ec 4.18
    mCH4bw = mCH4pw * np.exp(-2 * lambCH4 * m - (-3.89 * 10 ** (-3)) * m ** 2)

    # Ec. 4.23 Specific Volume of methane-free
    Vbw = 1 / (denbwnogas)

    # Step6: Calculate the formation volume factor

    # Ec 4.36

    if unit == 1:
        Bw = ((1000 + m * 58.44) * Vbw + mCH4bw * VMCH4B) / ((1000 + m * 58.44) * 0.99681786)
        Bw = Bw * 1

    else:
        Bw = ((1000 + m * 58.44) * Vbw + mCH4bw * VMCH4B) / ((1000 + m * 58.44) * 0.99681786)

    return Bw


def comp_bw_nogas(P, T, salinity, unit=1):
    """This function calculate the compressibility of methane-free brine
    The user must select UNITS 1 to Field units and 2 to International System

    IF UNITS 1 is selected introduce:
        | Pressure (P) in PSI
        | Temperature (T) in ºF
        | Salinity in ppm
        | Compressibility returns in psi-1

    IF UNITS 2 is selected introduce:
        | Pressure (P) in MPa
        | Temperature (T) in ºC
        | Salinity in ppm
        | Compressibility returns in MPa-1"""

    if unit == 1:
        P_array = np.array(P)
        P_array = P_array * (6.89475729 * 10 ** -3.0)
        T_array = np.array(T)
        T_array = (T_array - 32) * 5 / 9
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))
    else:
        P_array = np.array(P)
        T_array = np.array(T)
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    # Step 1: Calculate the compressibility coefficients

    # Ec 4.1 and coefficients Table 4_6
    a1_EW = 4.221
    a2_EW = -3.478
    a3_EW = 6.221
    a4_EW = 0.5182
    a5_EW = -0.4405

    a_T_EW = (a1_EW * (T_array / 100) ** 2 + a2_EW * (T_array / 100) + a3_EW) / (
                a4_EW * (T_array / 100) ** 2 + a5_EW * (T_array / 100) + 1)

    a1_FW = -11.403
    a2_FW = 29.932
    a3_FW = 27.952
    a4_FW = 0.20684
    a5_FW = 0.3768

    a_T_FW = (a1_FW * (T_array / 100) ** 2 + a2_FW * (T_array / 100) + a3_FW) / (
            a4_FW * (T_array / 100) ** 2 + a5_FW * (T_array / 100) + 1)

    # Ec 4.1 and coefficients Table 4_8

    a1_EM = 0
    a2_EM = 0
    a3_EM = 0.1249
    a4_EM = 0
    a5_EM = 0

    a_T_EM = (a1_EM * (T_array / 100) ** 2 + a2_EM * (T_array / 100) + a3_EM) / (
            a4_EM * (T_array / 100) ** 2 + a5_EM * (T_array / 100) + 1)

    a1_FM3_2 = -0.617
    a2_FM3_2 = -0.747
    a3_FM3_2 = -0.4339
    a4_FM3_2 = 0
    a5_FM3_2 = 10.26

    a_T_FM3_2 = (a1_FM3_2 * (T_array / 100) ** 2 + a2_FM3_2 * (T_array / 100) + a3_FM3_2) / (
            a4_FM3_2 * (T_array / 100) ** 2 + a5_FM3_2 * (T_array / 100) + 1)

    a1_FM1 = 0
    a2_FM1 = 9.917
    a3_FM1 = 5.1128
    a4_FM1 = 0
    a5_FM1 = 3.892

    a_T_FM1 = (a1_FM1 * (T_array / 100) ** 2 + a2_FM1 * (T_array / 100) + a3_FM1) / (
                a4_FM1 * (T_array / 100) ** 2 + a5_FM1 * (T_array / 100) + 1)

    a1_FM1_2 = 0.0365
    a2_FM1_2 = -0.0369
    a3_FM1_2 = 0
    a4_FM1_2 = 0
    a5_FM1_2 = 0

    a_T_FM1_2 = (a1_FM1_2 * (T_array / 100) ** 2 + a2_FM1_2 * (T_array / 100) + a3_FM1_2) / (
                a4_FM1_2 * (T_array / 100) ** 2 + a5_FM1_2 * (T_array / 100) + 1)

    # Ec 4.7
    Eb_T_M = a_T_EW + a_T_EM * m
    # Ec 4.8
    FB_T_M = a_T_FW + a_T_FM3_2 * m ** (3 / 2) + a_T_FM1 * m + a_T_FM1_2 * m ** (1 / 2)

    # Step 2: Calculate the Compressibility of methane-free Brine at (P,T,m)
    # Ec 4.9  (Compressibility of methane-free Brine at (P,T,m))

    if unit == 1:
        temp = 1 / 70 * (1 / (Eb_T_M * (P_array / 70) + FB_T_M))
        temp = temp * 6.89475729 * 10 ** (-3.0)

    else:
        temp = 1 / 70 * (1 / (Eb_T_M * (P_array / 70) + FB_T_M))

    return temp