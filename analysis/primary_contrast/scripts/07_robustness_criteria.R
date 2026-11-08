run_component("robustness_criteria", function() {
  need_primary()
  crit <- robustness_criteria(ST$primary$theta, ST$loo, ST$sub_primary$theta, ST$primary$R,
                              MONO, NK, CFG, ST$analysis$eight)
  crit$analysis <- ST$analysis$label
  write_result(crit, "robustness_summary.csv", RUN_ID)
})
