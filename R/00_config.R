# =====================================================================
# CNPS.cycle Shiny app - cycle configuration
# One row per sub-pathway, in the exact row order produced by
# *cyc.abundance(). Fields:
#   code    - short id (also the host-analysis function name in CNPS.cycle)
#   label   - exact Type label used in the abundance tables (and by the
#             fold.change heatmap `ylim()` subsetting trick)
#   cn      - Chinese display name for the UI
#   vx, vy  - viewport x/y for overlaying the mini heatmap on the cycle
#             background PDF (width = 0.035*Group_numb, height = 0.05)
#   host    - TRUE if the official script runs host-taxonomy analysis
#   title   - title drawn above the host.ratio plot
# =====================================================================

CNPS_PALETTE <- c("#B2182B", "#56B4E9", "#E69F00", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC",
                  "#66CC99", "#99999", "#ADD1E5")

CNPS_TAX_LEVELS <- c("门" = "Phylum", "纲" = "Class", "目" = "Order",
                     "科" = "Family", "属" = "Genus", "种" = "Species")

CNPS_CONFIG <- list(

  C = list(
    name = "碳循环", en = "Carbon", abun_fun = "Ccyc.abundance",
    host_fun = "Ccyc.host", bg_pdf = "Ccyc.pdf",
    heatmap_w = 2.3, heatmap_h = 3, fold_w = 2.3, fold_h = 4.5,
    processes = data.frame(
      code  = c("ACF", "ACH4O", "AR", "AnCF", "COo", "Fer", "Meth"),
      label = c("Aerobic C fixation", "Aerobic CH4 oxidation", "Aerobic respiration",
                "Anaerobic C fixation", "CO oxidation", "Fermentation", "Methanogenesis"),
      cn    = c("好氧固碳", "好氧甲烷氧化", "有氧呼吸", "厌氧固碳", "CO氧化", "发酵", "产甲烷"),
      vx    = c(0.692, 0.5935, 0.612, 0.322, 0.812, 0.36, 0.362),
      vy    = c(0.87, 0.23, 0.64, 0.86, 0.385, 0.45, 0.34),
      host  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
      title = c("Aerobic C fixation", "Aerobic CH4 oxidation", "Aerobic respiration",
                "Anaerobic C fixation", "CO oxidation", "Fermentation", "Methanogenesis"),
      stringsAsFactors = FALSE
    )
  ),

  N = list(
    name = "氮循环", en = "Nitrogen", abun_fun = "Ncyc.abundance",
    host_fun = "Ncyc.host", bg_pdf = "Ncyc.pdf",
    heatmap_w = 3.8, heatmap_h = 4, fold_w = 3.8, fold_h = 4.5,
    processes = data.frame(
      code = c("DNR.NO3.NO2", "DNR.NO2.NH4", "ANR.NO3.NO2", "ANR.NO2.NH4",
               "Denitrification.NO2.NO", "Denitrification.NO.N2O", "Denitrification.N2O.N2",
               "NF.N2.NH4", "Nitrification.NH4.NH2OH", "Nitrification.NH2OH.NO2",
               "Nitrification.NH2OH.NH4", "Anammox.NH4.N2H4", "Anammox.N2H4.N2",
               "NA.NH4.ON", "NM.ON.NH4", "NU.NH4", "NU.NO2", "NU.NO3"),
      label = c("Dissimilatory nitrate reduction, nitrate to nitrite",
                "Dissimilatory nitrate reduction, nitrite to ammonia",
                "Assimilatory nitrate reduction, nitrate to nitrite",
                "Assimilatory nitrate reduction, nitrite to ammonia",
                "Denitrification, nitrite to nitric oxide",
                "Denitrification, nitric oxide to nitrous oxide",
                "Denitrification, nitrous oxide to nitrogen",
                "Nitrogen fixation, nitrogen to ammonia",
                "Nitrification, ammonia to hydroxylamine",
                "Nitrification, hydroxylamine to nitrite",
                "Nitrification, hydroxylamine to ammonia",
                "Anammox, ammonia to hydrazine",
                "Anammox, hydrazine to nitrogen",
                "Nitrogen assimilation, ammonia to organ-N",
                "Nitrogen mineralization, organ-N to ammonia",
                "Ammonia uptake", "Nitrite uptake", "Nitrate uptake"),
      cn = c("异化硝酸盐还原 NO3->NO2", "异化硝酸盐还原 NO2->NH4",
             "同化硝酸盐还原 NO3->NO2", "同化硝酸盐还原 NO2->NH4",
             "反硝化 NO2->NO", "反硝化 NO->N2O", "反硝化 N2O->N2",
             "固氮 N2->NH4", "硝化 NH4->NH2OH", "硝化 NH2OH->NO2",
             "硝化 NH2OH->NH4", "厌氧氨氧化 NH4->N2H4", "厌氧氨氧化 N2H4->N2",
             "氮同化 NH4->有机氮", "氮矿化 有机氮->NH4",
             "氨吸收", "亚硝酸盐吸收", "硝酸盐吸收"),
      vx = c(0.383, 0.33, 0.583, 0.573, 0.165, 0.145, 0.183,
             0.583, 0.65, 0.648, 0.647, 0.423, 0.438,
             0.843, 0.753, 0.753, 0.763, 0.693),
      vy = c(0.085, 0.2, 0.12, 0.397, 0.3, 0.55, 0.84,
             0.9, 0.691, 0.347, 0.59, 0.55, 0.76,
             0.64, 0.87, 0.742, 0.28, 0.06),
      host = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE,
               TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, TRUE),
      title = c("Dissimilatory nitrate reduction (NO3- => NO2-)",
                "Dissimilatory nitrate reduction (NO2- => NH4+)",
                "Assimilatory nitrate reduction (NO3- => NO2-)",
                "Assimilatory nitrate reduction (NO2- => NH4+)",
                "Denitrification (NO2- => NO)",
                "Denitrification (NO => N2O)",
                "Denitrification (N2O => N2)",
                "Nitrogen fixation (N2 => NH4+)",
                "Nitrification (NH4+ => NH2OH)",
                "Nitrification (NH2OH => NO2-)",
                "Nitrification (NH2OH => NH4+)",
                "Anammox (NH4+ => N2H4)",
                "Anammox (N2H4 => N2)",
                "Nitrogen assimilation (NH4+ => organ-N)",
                "Nitrogen mineralization (organ-N => NH4+)",
                "Ammonia uptake", "Nitrite uptake", "Nitrate uptake"),
      stringsAsFactors = FALSE
    )
  ),

  P = list(
    name = "磷循环", en = "Phosphorus", abun_fun = "Pcyc.abundance",
    host_fun = "Pcyc.host", bg_pdf = "Pcyc.pdf",
    heatmap_w = 2.3, heatmap_h = 3, fold_w = 2.3, fold_h = 1.2,
    processes = data.frame(
      code  = c("RB", "X3"),
      label = c("PhoR-PhoB system", "SenX3-RegX3 system"),
      cn    = c("PhoR-PhoB 双组分系统", "SenX3-RegX3 双组分系统"),
      vx    = c(0.564, 0.408),
      vy    = c(0.725, 0.54),
      host  = c(TRUE, TRUE),
      title = c("PhoR-PhoB system", "SenX3-RegX3 system"),
      stringsAsFactors = FALSE
    )
  ),

  S = list(
    name = "硫循环", en = "Sulfur", abun_fun = "Scyc.abundance",
    host_fun = "Scyc.host", bg_pdf = "Scyc.pdf",
    heatmap_w = 3.2, heatmap_h = 3, fold_w = 3.2, fold_h = 4.5,
    processes = data.frame(
      code = c("ARS.SO4.APS", "ARS.APS.PAPS", "ARS.PAPS.SO3", "ARS.SO3.H2S",
               "DRS.SO4.APS", "DRS.APS.SO3", "DRS.SO3.H2S", "SOX",
               "SC.H2S.S", "SC.H2S.Sn", "SC.T.H2S", "SM",
               "SU.SO4", "SU.SO3", "SA"),
      label = c("Assimilatory sulfate reduction, sulfate to APS",
                "Assimilatory sulfate reduction, APS to PAPS",
                "Assimilatory sulfate reduction, PAPS to sulfite",
                "Assimilatory sulfate reduction, sulfite to sulfide",
                "Dissimilatory sulfate reduction and oxidation, sulfate to APS",
                "Dissimilatory sulfate reduction and oxidation, APS to sulfite",
                "Dissimilatory sulfate reduction and oxidation, sulfite to sulfide",
                "SOX system",
                "Sulfide cycling, sulfide to sulfur",
                "Sulfide cycling, sulfide to (sulfide)n",
                "Sulfide cycling, thisulfate to sulfide",
                "Sulfur mineralization",
                "Sulfate uptake", "Sulfite uptake", "Sulfur assimilation"),
      cn = c("同化硫酸盐还原 SO4->APS", "同化硫酸盐还原 APS->PAPS",
             "同化硫酸盐还原 PAPS->SO3", "同化硫酸盐还原 SO3->H2S",
             "异化硫酸盐还原氧化 SO4->APS", "异化硫酸盐还原氧化 APS->SO3",
             "异化硫酸盐还原氧化 SO3->H2S", "SOX 系统",
             "硫循环 H2S->S", "硫循环 H2S->(H2S)n", "硫循环 硫代硫酸盐->H2S",
             "硫矿化", "硫酸盐吸收", "亚硫酸盐吸收", "硫同化"),
      vx = c(0.225, 0.56, 0.69, 0.725, 0.223, 0.475, 0.73,
             0.335, 0.73, 0.59, 0.61, 0.872, 0.18, 0.81, 0.47),
      vy = c(0.6, 0.07, 0.25, 0.55, 0.13, 0.15, 0.66,
             0.89, 0.845, 0.79, 0.895, 0.73, 0.83, 0.34, 0.55),
      host = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
               TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
      title = c("Assimilatory sulfate reduction (SO4 => APS)",
                "Assimilatory sulfate reduction (APS => PAPS)",
                "Assimilatory sulfate reduction (PAPS => SO3)",
                "Assimilatory sulfate reduction (SO3 => H2S)",
                "Dissimilatory sulfate reduction (SO4 => APS)",
                "Dissimilatory sulfate reduction (APS => SO3)",
                "Dissimilatory sulfate reduction (SO3 => H2S)",
                "SOX system",
                "Sulfide cycling (H2S => S)",
                "Sulfide cycling (H2S => (H2S)n)",
                "Sulfide cycling (thiosulfate => H2S)",
                "Sulfur mineralization",
                "Sulfate uptake", "Sulfite uptake", "Sulfur assimilation"),
      stringsAsFactors = FALSE
    )
  )
)

CNPS_CYCLE_KEYS <- c("C", "N", "P", "S")

# KO lists per cycle (official script, chunks 4 / 11 / 18 / 25)
CNPS_CYCLE_KOS <- list(
  C = c("K00855","K01061","K01602","K08684","K02256","K02262","K02274","K02276",
        "K00174","K00175","K00244","K01648","K00194","K00197","K03518","K03519",
        "K03520","K00016","K00400","K00401"),
  N = c("K00370","K00371","K00374","K02567","K02568","K00362","K00363","K03385",
        "K15876","K00367","K00372","K00366","K00368","K15864","K04561","K02305",
        "K00376","K02588","K02586","K02591","K10944","K10945","K10946","K10535",
        "K05601","K20932","K20933","K20934","K20935","K01915","K00265","K00266",
        "K00264","K00284","K00260","K15371","K00261","K00262","K01455","K01501",
        "K01725","K00926","K00549","K02575","K15576","K15577","K15578","K15579"),
  P = c("K07636","K07657","K07658","K07768","K07776"),
  S = c("K00956","K00957","K00955","K00860","K00390","K00380","K00381","K00392",
        "K00958","K00988","K00394","K00395","K11180","K11181","K17222","K17223",
        "K17224","K17225","K17226","K17227","K17229","K17230","K17218","K08352",
        "K08354","K04091","K00299","K16968","K16969","K15762","K15765","K03119",
        "K00456","K17217","K02045","K02046","K02047","K02048","K15551","K15552",
        "K10831","K15553","K15554","K15555","K01739","K10764","K01738","K17228")
)
