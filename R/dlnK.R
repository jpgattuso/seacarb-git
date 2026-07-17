# This file is part of seacarb.
#
# dlnK: analytic logarithmic derivatives  dln(K)/dT  and  dln(K)/dS  of the
#       equilibrium constants used by buffderiv().
#
# seacarb supplies every K VALUE (via K1(), K2(), Kb(), Kw(), ... ), so nothing
# here recomputes a constant.  What seacarb does not supply is dK/dT and dK/dS,
# and that is all this file provides.  The caller reassembles them as
#
#   total-scale-native (K1, K2, Kb):   dK/dY = K * dlnK/dY
#   SWS-native  (Kw, Ksi, K1p..K3p):   dK/dY = K * ( dlnKsws/dY + dlnSig/dY )
#
# where Sig = kSWS2total.  Ks and Kf enter only through Sig.
#
# These are the exact derivatives of the formulations selected by
# k1k2="l", ks="d", kf="dg", b="u74", pHscale="T"  (Orr et al. 2015).
# They are validated in validate_buffderiv.R against finite differences of
# seacarb's OWN K functions, so they cannot silently drift out of sync.

dlnK <- function(S, T, Ks, Kf) {
    TK <- T + 273.15
    lnTK <- log(TK); sqS <- sqrt(S); ln10 <- log(10)
    Io <- 19.924*S/(1000 - 1.005*S); sqIo <- sqrt(Io)
    Cl <- S/1.80655
    ST <- 0.14*Cl/96.062
    FT <- 6.7e-5*Cl/18.9984
    dIo_dS  <- 19.924*1000/(1000 - 1.005*S)^2
    dlgS_dS <- -0.001005/(1 - 0.001005*S)   # d/dS of log(1 - 0.001005*S)

    ## --- K1, K2: Lueker et al. (2000), log10 form, total scale --------------
    K1_T <- ln10*( 3633.86/TK^2 - 9.67770/TK)
    K1_S <- ln10*( 0.011555 - 2*0.0001152*S)
    K2_T <- ln10*( 471.78/TK^2 + 3.16967/TK)
    K2_S <- ln10*( 0.01781 - 2*0.0001122*S)

    ## --- Kb: Dickson (1990b), total scale ------------------------------------
    tb1     <- -8966.90 - 2890.53*sqS - 77.942*S + 1.728*S^1.5 - 0.0996*S*S
    dtb1_dS <- -2890.53*0.5/sqS - 77.942 + 1.728*1.5*sqS - 0.1992*S
    Kb_T <- -tb1/TK^2 + (-24.4344 - 25.085*sqS - 0.2474*S)/TK + 0.053105*sqS
    Kb_S <- dtb1_dS/TK + 137.1942*0.5/sqS + 1.62142 +
            (-25.085*0.5/sqS - 0.2474)*lnTK + 0.053105*0.5/sqS*TK

    ## --- Ks: Dickson (1990b), free scale -------------------------------------
    A_ <- -13856/TK + 324.57 - 47.986*lnTK
    B_ <-  35474/TK - 771.54 + 114.723*lnTK
    Ks_T <- 4276.1/TK^2 - 23.093/TK + (13856/TK^2 - 47.986/TK)*sqIo +
            (-35474/TK^2 + 114.723/TK)*Io +
            (2698/TK^2)*Io*sqIo - (1776/TK^2)*Io*Io
    Ks_S <- (A_*0.5/sqIo + B_ + (-2698/TK)*1.5*sqIo + (1776/TK)*2*Io)*dIo_dS + dlgS_dS

    ## --- Kf: Dickson & Riley (1979), free scale ------------------------------
    Kf_T <- -1590.2/TK^2
    Kf_S <- 1.525*0.5/sqIo*dIo_dS + dlgS_dS

    ## --- Sig = kSWS2total = (1 + ST/Ks) / (1 + ST/Ks + FT/Kf) ----------------
    U <- 1 + ST/Ks
    V <- U + FT/Kf
    dST_dS <- 0.14/96.062/1.80655
    dFT_dS <- 6.7e-5/18.9984/1.80655
    dU_dT <- -(ST/Ks)*Ks_T
    dV_dT <- dU_dT - (FT/Kf)*Kf_T
    dU_dS <- dST_dS/Ks - (ST/Ks)*Ks_S
    dV_dS <- dU_dS + dFT_dS/Kf - (FT/Kf)*Kf_S
    Sig_T <- dU_dT/U - dV_dT/V
    Sig_S <- dU_dS/U - dV_dS/V

    ## --- Sig_x: the SWS->total factor that seacarb's Kn() ACTUALLY uses -------
    # Kn() and Khs() do NOT accept a kSWS2scale argument (unlike Kw, Ksi, K1p...),
    # so they call kconv() with its default kf="x".  Kf.R then does
    #     is_outrange <- T>33 | T<10 | S<10 | S>40
    #     kf[is_x] <- "pf" ; kf[is_x & is_outrange] <- "dg"
    # i.e. Perez & Fraga IN range, Dickson & Riley outside.  Khs is total-scale
    # native so the factor round-trips and cancels, but Kn is SWS native and does
    # NOT cancel.  To differentiate seacarb's Kn() we must follow the active branch.
    #
    # WARNING: this makes dKn/dT and dKn/dS DISCONTINUOUS across T=10, T=33, S=10,
    # S=40.  Whenever NH4t != 0, dBe/dT and dBe/dS inherit that discontinuity.  It is
    # a seacarb API defect, not a property of the chemistry.  The clean fix upstream
    # is to give Kn() and Khs() a kSWS2scale argument like the other constants have.
    in_pf <- !(T > 33 | T < 10 | S < 10 | S > 40)
    daph_dT <- -(ST/Ks)*Ks_T
    daph_dS <- dST_dS/Ks - (ST/Ks)*Ks_S
    Kfpf   <- exp(874/TK - 9.68 + 0.111*sqrt(S))/U        # free scale; U = 1 + ST/Ks
    Kfpf_T <- -874/TK^2        - daph_dT/U
    Kfpf_S <- 0.111*0.5/sqrt(S) - daph_dS/U
    Kfx   <- ifelse(in_pf, Kfpf,   Kf)
    Kfx_T <- ifelse(in_pf, Kfpf_T, Kf_T)
    Kfx_S <- ifelse(in_pf, Kfpf_S, Kf_S)
    Vx <- U + FT/Kfx
    dVx_dT <- dU_dT - (FT/Kfx)*Kfx_T
    dVx_dS <- dU_dS + dFT_dS/Kfx - (FT/Kfx)*Kfx_S
    Sigx_T <- dU_dT/U - dVx_dT/Vx
    Sigx_S <- dU_dS/U - dVx_dS/Vx

    ## --- SWS-native constants: Millero (1995) --------------------------------
    Kw_T <- 13847.26/TK^2 - 23.6521/TK + (-118.67/TK^2 + 1.0495/TK)*sqS
    Kw_S <- (118.67/TK - 5.977 + 1.0495*lnTK)*0.5/sqS - 0.01615

    Bsi   <- -8904.2 - 458.79*sqIo + 188.74*Io - 12.1652*Io*Io
    Ksi_T <- -Bsi/TK^2 - 19.334/TK
    Ksi_S <- (3.5913*0.5/sqIo - 1.5998 + 2*0.07871*Io +
              (-458.79*0.5/sqIo + 188.74 - 2*12.1652*Io)/TK)*dIo_dS + dlgS_dS

    K1p_T <- 4576.752/TK^2 - 18.453/TK + (106.736/TK^2)*sqS + (0.65643/TK^2)*S
    K1p_S <- (-106.736/TK + 0.69171)*0.5/sqS + (-0.65643/TK - 0.01844)
    K2p_T <- 8814.715/TK^2 - 27.927/TK + (160.34/TK^2)*sqS + (-0.37335/TK^2)*S
    K2p_S <- (-160.34/TK + 1.3566)*0.5/sqS + (0.37335/TK - 0.05778)
    K3p_T <- 3070.75/TK^2 + (-17.27039/TK^2)*sqS + (44.99486/TK^2)*S
    K3p_S <- (17.27039/TK + 2.81197)*0.5/sqS + (-44.99486/TK - 0.09984)

    ## --- K2si: Millero (1995), log10, TOTAL scale native --------------------
    # 10^logK * ktotal2SWS_P0 * kSWS2total = 10^logK exactly at P = 0.
    # logK depends on T only, so dK2si/dS = 0 on the total scale.
    K2si_T <- ln10*( 4465.18/TK^2 - 0.021952 )
    K2si_S <- 0

    ## --- Khs: TOTAL scale native (same round trip as Kb, K2si) ---------------
    Khs_T <- 13275.3/TK^2 - 34.6435/TK
    Khs_S <- 0.3449*0.5/sqrt(S) - 0.0274

    ## --- Kn: SWS native, so the caller must add dlnSig/dY ---------------------
    Kn_T <- 6285.33/TK^2 + 0.0001635 + (123.7184/TK^2)*sqrt(S) - (3.17556/TK^2)*S
    Kn_S <- (0.46532 - 123.7184/TK)*0.5/sqrt(S) + (-0.01992 + 3.17556/TK)

    list(Sigx_T=Sigx_T, Sigx_S=Sigx_S,
         K2si_T=K2si_T, K2si_S=K2si_S, Khs_T=Khs_T, Khs_S=Khs_S,
         Kn_T=Kn_T, Kn_S=Kn_S,
         K1_T=K1_T, K1_S=K1_S, K2_T=K2_T, K2_S=K2_S, Kb_T=Kb_T, Kb_S=Kb_S,
         Ks_T=Ks_T, Ks_S=Ks_S, Kf_T=Kf_T, Kf_S=Kf_S, Sig_T=Sig_T, Sig_S=Sig_S,
         Kw_T=Kw_T, Kw_S=Kw_S, Ksi_T=Ksi_T, Ksi_S=Ksi_S,
         K1p_T=K1p_T, K1p_S=K1p_S, K2p_T=K2p_T, K2p_S=K2p_S,
         K3p_T=K3p_T, K3p_S=K3p_S)
}
