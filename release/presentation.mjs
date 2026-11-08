// Render the seven-slide research talk with @oai/artifact-tool or Artifact Session.
export function buildPresentation(presentation) {
  const font = 'Arial';
  presentation.theme.defaultFont = font;
  function text(slide, value, x, y, w, h, size = 26, color = '#25394A', bold = false) {
    const shape = slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
    shape.text = value;
    shape.text.style = {typeface:font,fontSize:size,color,bold};
    return shape;
  }
  function slide(title, notes) {
    const s = presentation.slides.add({width:1280,height:720});
    s.background.fill = '#FAFAF7';
    text(s,title,64,42,1152,100,40,'#142C3D',true);
    s.speakerNotes.text = notes;
    return s;
  }
  function table(s, values, y, height, widths) {
    const t = s.tables.add({rows:values.length,columns:values[0].length,left:64,top:y,width:1152,height,columnWidths:widths,values});
    t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({textStyle:{typeface:font,fontSize:24,color:'#25394A'},fill:'#FAFAF7',margins:{left:12,right:12,top:12,bottom:12}});
    t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({textStyle:{typeface:font,fontSize:24,color:'#FFFFFF',bold:true},fill:'#205A69'});
    t.borders.assign({fill:'#CED9DC',width:1,style:'solid'});
    return t;
  }
  let s=slide('Context dependence and measurement sensitivity', 'Full study title: Context Dependence and Measurement Sensitivity of Interferon-β Induction Contrasts in Human Monocytes and NK Cells: A Three-Dataset Computational Study. Zeyu Yao. v1.2.0. No new experimental data were generated. Sources: release/research_report.md; README.md.');
  text(s,'Interferon-β induction contrasts in human monocytes and NK cells',64,184,1110,140,44,'#205A69',true);
  text(s,'A three-dataset computational study',64,355,1100,58,30);
  text(s,'Transferability of a fixed IFN-response contrast requires evaluation across experimental settings.',64,452,1050,110,34);
  text(s,'Zeyu Yao  ·  8 November 2026  ·  v1.2.0',64,642,1110,36,22);

  s=slide('Biological motivation and statistical question','Panel from Rigby et al., https://doi.org/10.1038/s44319-026-00750-3. Dataset sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226572 . Counts summed by donor, lineage and condition; separate TMM normalization by lineage; mean ten-gene stimulated-minus-control log2 CPM; equal weight per donor. Primary prior count 2. See README Study design.');
  text(s,'IFN genes can be induced in both cell types while their induction contrast depends on measurement and experimental setting.',64,148,1130,90,28);
  text(s,'δdonor = monocyte induction − NK induction       θ = mean donor δ',64,260,1152,55,28,'#205A69',true);
  table(s,[['Dataset','Donors / preparation','IFN-β exposure'],['GSE96583','8 SLE / mixed PBMC','100 U/mL · 6 h'],['HIRISA / GSE306664','5 healthy / enriched cells','100 U/mL · 21 h'],['GSE226572','3 healthy / mixed PBMC','1,000 U/mL · 4–24 h']],346,246,[330,480,342]);
  text(s,'Donors are the replicates. These are relative transcript measurements.',64,628,1152,52,25);

  s=slide('GSE96583: measurement sensitivity','Sources: results/primary_contrast/02_primary_estimate/primary.csv; results/primary_contrast/14_fixed_offset_measurement/fixed_offset_summary.csv; results/primary_contrast/08_pseudocount_sensitivity/pc_library_sizes.csv. Conditional primary 95% t interval 1.163–2.506. Sensitivities reuse the same donors; changing the CPM offset redistributes regularization with the preferred measurement specification unresolved.');
  table(s,[['Expression measurement','Mean contrast θ','Positive donors'],['edgeR prior count 2','1.834','8 / 8'],['edgeR prior count 0.5','1.595','7 / 8'],['log₂(CPM + 1)','1.233','7 / 8']],160,286,[600,275,277]);
  text(s,'The fixed offset lowers θ by 0.601 and reverses one donor’s sign.',64,484,1140,65,29,'#205A69',true);
  text(s,'Prior count 2 adds 17.01 CPM in smaller NK libraries versus 2.52 CPM in monocytes. Sparse control NK genes drive much of the sensitivity.',64,574,1140,95,27);

  s=slide('HIRISA: the contrast reverses across context','Sources: results/independent_cohort/h2_theta.csv; h2_offset_sensitivity.csv; results/independent_cohort/triangulation/core_gene_comparison.csv and summary_metrics.csv. HIRISA source https://doi.org/10.21203/rs.3.rs-8197575/v1. MAST summaries are medians across donor-level fits on a different scale from θ. Both analyses use the same cells and assess analytical robustness within one cohort. Broader Hallmark set: mean −0.990, four of five donors negative, interval includes zero.');
  table(s,[['Ten-gene panel','GSE96583','HIRISA'],['Mean contrast θ','+1.834','−4.609'],['Negative donors','0 / 8','5 / 5'],['Change with fixed offset','−0.601','+0.155']],160,276,[600,275,277]);
  text(s,'Author MAST summaries agree: greater NK induction for all 10 genes.',64,470,1140,76,29,'#205A69',true);
  text(s,'High control monocyte expression coincides with weak induction in 4/5 donors. Health, preparation, duration and assay differ together; cause remains unresolved.',64,573,1140,104,27);

  s=slide('GSE226572: the contrast rises from 4 h to 24 h','Source: results/time_course/g3_gamma_by_donor.csv and results/time_course/g3_sensitivity.csv, primary marker annotation, edgeR prior count 2. Mean θ4h −0.173; θ24h 0.382; Γ 0.555. Three donors; smallest two-sided sign-flip p = 0.25. Common untreated reference cancels in the temporal contrast algebraically, while normalization can depend on reference profile composition. Original response dynamics: de Cevins et al., https://pmc.ncbi.nlm.nih.gov/articles/PMC10772457/ Figure 2.');
  const chart = s.charts.add('line',{position:{left:70,top:162,width:690,height:410},categories:['4 h','24 h'],series:[{name:'D34',values:[-0.3858653753,0.2842731182],fill:'#205A69'},{name:'D38',values:[0.3033194959,0.7450526976],fill:'#B87532'},{name:'D39',values:[-0.4375845205,0.1159100768],fill:'#8071A0'}],hasLegend:true,yAxis:{title:'Monocyte − NK contrast (log₂ units)',minimumScale:-0.6,maximumScale:0.9,majorGridlines:{fill:'#D8DEDF',width:1}},xAxis:{title:'IFN-β exposure'}});
  chart.xAxis = {title:'',tickLabelPosition:'low',textStyle:{fontSize:20}};
  chart.yAxis = {title:{text:'Monocyte − NK contrast (log₂ units)',textStyle:{fontSize:20}},textStyle:{fontSize:18}};
  chart.legend = {position:'bottom',textStyle:{fontSize:20}};
  for (const style of [chart.legend.textStyle, chart.dataLabels.textStyle, chart.xAxis.textStyle, chart.yAxis.textStyle, chart.yAxis.title.textStyle]) style.typeface = font;
  text(s,'Γ = +0.555',810,192,370,75,42,'#205A69',true);
  text(s,'Positive change in all three donors',810,294,370,95,30);
  text(s,'The observed direction opposed the prespecified decrease.',810,434,370,112,29);
  text(s,'A simple duration-only explanation is challenged; duration × culture interactions remain possible.',64,609,1152,70,27);

  s=slide('The temporal increase survives four annotations','Source: results/time_course/annotation_validation/gamma_by_annotation.csv, prior count 2 rows. All four mean changes and positive-donor counts are read from saved results. Original and expanded marker rules; independent depositor classical and broader monocyte labels. Metadata mapped for all 115,503 author-labelled cells. The four definitions assess annotation sensitivity within the same three donors.');
  table(s,[['Monocyte definition','Mean change Γ','Positive donors'],['Marker CD14','0.555','3 / 3'],['Marker CD14 / CD16','0.595','3 / 3'],['Depositor classical','0.597','3 / 3'],['Depositor broader','0.553','3 / 3']],160,350,[600,275,277]);
  text(s,'Independent labels support annotation robustness within the same three donors.',64,549,1140,72,29,'#205A69',true);
  text(s,'Contribution: a fixed donor-level contrast and its robustness. IFN-response dynamics were already studied in the original publication.',64,628,1140,60,23);

  s=slide('Methodological lessons and limits','Sources: README Conclusions and Limitations; release/FINALIZATION.md; release/research_report.md. No causal role established for enrichment, baseline activation, disease status or a single mechanism. GSE96583 processing pool is confounded with stimulation. GSE226572 second untreated library role undocumented. Historical dirty run IDs retained. No expensive analysis rerun during finalization. Future spatial-omics investigation is a proposed distinct project and has not been executed.');
  text(s,'Report direction, magnitude and donor-level stability separately.',64,160,1150,70,33,'#205A69',true);
  text(s,'Small donor counts (8 / 5 / 3), different designs and unmeasured mechanisms limit generalization. Automated checks assess internal consistency; scientific validity depends on study design and interpretation.',64,268,1130,130,29);
  text(s,'Completed artifact: locked specifications, scripts, figures, tables, checksums, software versions and reproduction instructions.',64,430,1130,103,28);
  text(s,'Close this investigation. A distinct spatial-omics or cell-neighbourhood robustness question is a possible next project.',64,574,1130,108,28);
  if (presentation.slides.items.length !== 7) throw new Error('Expected seven slides');
}
