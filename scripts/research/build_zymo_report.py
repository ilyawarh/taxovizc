"""Build the Zymo publication and standalone aggregate context; no inference rerun."""
from pathlib import Path
import csv,json,re,html,hashlib
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from PIL import Image as PILImage
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
R=Path(__file__).resolve().parents[2]; D=R/'runs/zymo_fecal_20261005/derived'; E=R/'docs/computing/evidence/zymo_20261005'; F=R/'docs/research/figures/zymo_20261005'
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def js(p):return json.loads(p.read_text(encoding='utf-8'))
review=(R/'docs/research/ZYMO_FECAL_REVIEW_20261005.md').read_text(encoding='utf-8')
experiment=(R/'docs/computing/experiments/EXP-20261005-ZYMO.md').read_text(encoding='utf-8')
supp=js(E/'supplemental_findings.json'); bins=rows(E/'bin_summary.csv')
extra='''## Additional findings and interpretation

Length is associated with classifier agreement: K=C species resolution rises from 52.6% for reads up to 2 kb to 97.1% above 30 kb. This is not a measured error rate or a causal length effect; taxon and length distributions are confounded. The independent length and GC panels help distinguish these associations from taxonomic islands.

K=C species labels resolve 199,035 reads (68.58%); 91,192 reads remain unresolved (31.42% of reads, 15.17% of bases). B. fragilis and E. faecalis together account for 35.67% of reads but 61.29% of bases. Read counts and base counts therefore produce substantially different pictures of dominance. Neither equals cell abundance.

P. gingivalis is detected sparsely: 42 MetaPhlAn, 407 KrakenUniq and 406 Centrifuger calls. Gemella-related genus calls are also sparse (22, 2 and 19 respectively); Centrifuger has only two G. morbillorum species calls. These are per-tool observations, not confirmed counts of true organisms, and sparse calls do not establish absence relative to the workbook.

Bin 0 illustrates why a coherent island cannot resolve conflicting taxonomy. KrakenUniq labels 13,368 of 13,611 reads P. anaerobius (98.2%), MetaPhlAn 3,757, and Centrifuger only 338. Centrifuger additionally assigns P. equinus (1,262), Caudoviricetes (1,078), P. porci (955), or unclassified (3,276). Independent sequence evidence is needed before reassignment.

Bin 7 contains 16,557 reads: Centrifuger assigns E. coli to 8,674 and lambda phage to 4,792, whereas KrakenUniq assigns E. coli to 8,105 and MetaPhlAn leaves 16,553 unclassified at species level. Related ambiguous bins 11 and 13 contain 981 and 1,028 reads. Bin 11 includes 557 Centrifuger lambda calls; bin 13 includes 678 E. coli calls. Mobile sequence, controls, reference overlap and parser effects remain competing hypotheses, not demonstrated contamination.

LRBinner's earlier attempts yielded 8, 15 and 19 bins before the final 23-bin result. The final resumed invocation took 52,249.43 seconds (14.5 hours), excluding preprocessing and training. Exhaustive cluster searching followed by assignment of remaining reads explains complete bin coverage; coverage alone does not measure bin purity.

## qPCR audit detail

The pure-Zymo block reports universal 16S Cq 19.87 / 1,400,000 copies, BFT 39.88 / 1.59, Porphyr 39.04 / 2.8, PKS 38.12 / 7.98 and Gem 36.67 / 25.1. The separate culture-associated block reports 16S 18.60 / 3,302,784, BFT 20.64 / 716,909, ENTRC 20.09 / 1,393,727, Fus 24.44 / 89,731, Gem 26.55 / 16,500, Pep 23.07 / 286,828, PKS 31.87 / 578, Porphyr 21.72 / 256,648 and Strep 24.27 / 20,188. These are workbook values; the second block is not explicitly identified as the sequenced mixture.

The common Sheet 2 formula combines H*10/5 + 6*O*1000/5, but ENTRC and Fus omit the factor 6. Preparation notes do not establish which convention is correct. All 39 target formulas, cached results, descriptions and cell locations are preserved in the JSON handoff. They were read, not recalculated or silently corrected. No toxin carriage, strain identity or absolute organism concentration is validated.
'''
context={'schema':'smrl.zymo.agent-context.v1','prepared':'2026-10-06','scope':'Completed descriptive review; no model rerun, taxonomy correction or independent accuracy validation','user_clarifications':['Intentional custom ontime-based earliest-sequence 3% selection; not random sampling and not a sampling bug.','Product, lot, dilution and additions cannot be clarified beyond the supplied workbook.','Current notebook was saved and fetched after initial stale snapshot.'],'report_markdown':review,'additional_analysis_markdown':extra,'experiment_record_markdown':experiment,'evidence':{p.stem:js(p) for p in E.glob('*.json') if p.name!='publication_manifest.json'},'all_bins':bins,'qpcr_all_targets':rows(D/'qpcr_formulas.csv'),'workbook_16S_reference':rows(D/'workbook_reference_16S.csv'),'taxon_profiles':{p.stem:rows(p) for p in D.glob('taxa_*.csv')},'limitations':['No independently validated read-level truth','Complete v1 fitted PCA state not retained','Classifier databases/taxonomy hashes and historical LRBinner seeds unverified','Full FASTQ checksum not computed; sequence-quality effects not assessed','Sampler timing fields and code not audited; file position is not acquisition time','No strain, toxin or contamination confirmation'],'artifact_policy':'Large originals, matrices and per-read IDs stay in ignored runs/zymo_fecal_20261005; this context contains aggregate evidence only.'}
cp=R/'docs/research/ZYMO_AGENT_CONTEXT_20261006.json';cp.write_text(json.dumps(context,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
for name,file in [('Arial','arial.ttf'),('Arial-Bold','arialbd.ttf'),('Arial-Italic','ariali.ttf')]:pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='Arial-Bold',italic='Arial-Italic',boldItalic='Arial-Bold')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Body',fontName='Arial',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#233343')))
styles.add(ParagraphStyle(name='SmallText',parent=styles['Body'],fontSize=8,leading=10,spaceAfter=3))
for n,size in [('Title',26),('Heading1',19),('Heading2',14)]:styles[n].fontName='Arial-Bold';styles[n].fontSize=size;styles[n].leading=size+5;styles[n].textColor=colors.HexColor('#124D65')
story=[]
def fmt(s):
    s=s.replace('в†’',' to ').replace('вЂ“','-').replace('вЂ”','-').replace('вЂ‘','-')
    s=html.escape(s)
    s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',lambda m:('<link href="'+m[2]+'" color="#167694">'+m[1]+'</link>') if m[2].startswith('http') else m[1],s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s);s=re.sub(r'(?<!\*)\*([^*]+)\*',r'<i>\1</i>',s)
    return s.replace('`','')
def p(s,style='Body'):return Paragraph(fmt(s),styles[style])
def table(data,widths=None):
    t=Table([[p(str(c),'SmallText') for c in row] for row in data],colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DBEAF0')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#6494A7')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F4F7F8')])]))
    story.extend([t,Spacer(1,10)])
def md(text):
    lines=text.splitlines();i=0
    while i<len(lines):
        line=lines[i].strip();i+=1
        if not line:continue
        if line.startswith('|'):
            data=[]
            while True:
                if not re.match(r'^\|[\s:|\-]+\|$',line):data.append([c.strip() for c in line.strip('|').split('|')])
                if i>=len(lines) or not lines[i].strip().startswith('|'):break
                line=lines[i].strip();i+=1
            table(data,[45,85,75,294] if data[0][0]=='LRBinner bin' else [499/len(data[0])]*len(data[0]));continue
        if line.startswith('# '):story.append(p(line[2:],'Heading1'))
        elif line.startswith('## '):story.append(p(line[3:],'Heading2'))
        else:story.append(p(line))
def page(title):story.extend([PageBreak(),p(title,'Heading1')])
story.extend([p('Zymo fecal experiment','Title'),p('TaxoViz maps, LRBinner bins and qPCR evidence','Heading2'),p('SMRL research report | 6 October 2026'),Spacer(1,14),p('The dominant groups are consistent with several added-culture candidates. The strongest unresolved findings are classifier conflict in bin 0, E. coli/lambda-related ambiguity in bins 7/11/13, and the association of map shape with read length. Version 2 changes the display; it has not demonstrated improved classification.'),p('Scope: saved-run analysis of 290,227 intentionally selected early sequences, linked to all 9,674,223 reads in the LRBinner run. No labels were corrected and no embedding or model was rerun.'),p('How to read the evidence','Heading2'),p('Observed: file checks, read joins, plotted coordinates and saved labels. Stored: notebook settings and historical logs. Hypotheses: biological and representation-based explanations for structure. Proposed: independent validation and correctness fixes. Tool agreement and bin membership are diagnostic evidence, not ground truth.'),p('Contents','Heading2'),p('Research interpretation and reference comparison; original-coordinate figures for bins, taxa, length and GC; quantitative neighborhood diagnostics; additional findings and qPCR audit; complete bin inventory; methods, provenance and sources.'),p('Companion: ZYMO_AGENT_CONTEXT_20261006.json contains the report text, complete aggregate evidence, all 23 bins, all 39 qPCR target records and taxonomic profiles. Large originals and read identifiers remain in the local ignored run archive.')])
page('Research interpretation');md(review.split('\n',1)[1])
for title,file,caption in [('LRBinner bins on the TaxoViz maps','bin_comparison.png','All 290,227 reads use the original HTML coordinates. Colors identify the 23 LRBinner bins. Small bins remain drawn even where text labels are suppressed. Several bins share broad regions, and bin membership does not establish species identity.'),('Taxonomic agreement on the same coordinates','taxa_comparison.png','K=C denotes identical KrakenUniq and Centrifuger species strings, excluding unclassified. The leading eight species are colored; other resolved and unresolved calls remain visible. Numbers refer to the legend, not LRBinner bin IDs.'),('Read length and the petal structures','length_comparison.png','Shared log10 length scale across both methods and versions. Many peripheral petals contain longer reads, whereas dense cores often contain shorter reads. Shape alone cannot identify strains, plasmids or error classes.'),('GC composition across versions','gc_comparison.png','Shared GC scale across both versions and methods. Major low- and high-GC islands recur, while geometry changes. GC is part of the TaxoViz feature input, so its organization is not independent validation.')]:
    page(title);iw,ih=PILImage.open(F/file).size;im=Image(str(F/file),width=499,height=499*ih/iw);story.extend([im,Spacer(1,12),p(caption)])
page('Neighborhood diagnostics')
pairs=[['Map','Same bin','Same K=C species*','Median |log2 length ratio|','Median GC difference'],['t-SNE v1','88.98%','90.46%','0.2497','0.00566'],['t-SNE v2','89.15%','90.02%','0.2342','0.00555'],['UMAP v1','85.54%','85.82%','0.3510','0.00928'],['UMAP v2','86.43%','86.39%','0.2886','0.00818']]
table(pairs,[85,80,110,120,104]);story.append(p('10,000 deterministic query reads, 15 exact Euclidean neighbors in each saved 2D map, seed 42; self excluded. *Only pairs resolved at both ends contribute to species agreement: 90,723 / 90,448 / 88,183 / 88,480 pairs respectively. Other diagnostics use all 150,000 directed pairs.'))
story.append(p('Exact cross-version neighbor overlap is 47.13% for t-SNE and 11.65% for UMAP. Cross-method overlap is 15.99% in v1 and 19.20% in v2. Broad concordance can coexist with substantially changed fine neighborhoods. No confidence interval, accuracy claim or original-space preservation test is implied.'))
table([['Length interval (bp)','Reads','K=C species resolved']]+[[f"{x['lower_exclusive_bp']:,} < L <= {x['upper_inclusive_bp']:,}",f"{x['reads']:,}",f"{100*x['KC_species_resolved_fraction']:.2f}%"] for x in supp['length_strata']],[230,120,149])
story.append(p('Legacy agreement categories hide conflicts: 37,289 of 43,794 species K_only reads and 5,617 of 6,087 M_only reads contain conflicting classified calls. Even 29,310 of 166,911 KC-category species reads have a conflicting third-tool label. Three-tool agreement is 20.94% at genus and 11.07% at species; it is not calibrated confidence.'))
page('Additional findings and qPCR');md(extra)
page('Complete 23-bin inventory')
story.append(p('Ordered by plotted size. Taxon percentage uses every plotted read in the bin, including unresolved reads. Selection % is plotted/full-bin reads. Top taxon is K=C species agreement, not a truth or purity estimate.'))
table([['Bin','Full reads','Plotted','Selection %','Top K=C species','% of bin']]+[[b['bin'],f"{int(b['full_reads']):,}",f"{int(b['plot_reads']):,}",f"{float(b['sampling_pct']):.2f}",b['top_KC_species'],f"{float(b['top_species_pct_all_bin_reads']):.1f}"] for b in bins],[30,75,60,62,215,57])
page('Methods, provenance and sources');md(experiment.split('\n',1)[1].replace('No commit, merge, push, production-code edit or remote training was performed.','At the close of the October 5 analysis, no commit, merge, push, production-code edit or remote training had been performed. Report publication and Git synchronization are a separate October 6 delivery step.').replace('Two diagnostic taxon labels can overlap at nearby medians; CSVs preserve exact identities.','The October 6 publication separates the nearby taxon annotations with leader lines; CSVs preserve exact identities.'))
out=R/'docs/research/reports/Zymo_Experiment_20261006.pdf';out.parent.mkdir(parents=True,exist_ok=True)
def footer(c,d):
    c.setStrokeColor(colors.HexColor('#CADCE4'));c.line(48,43,547,43);c.setFont('Arial',8);c.setFillColor(colors.HexColor('#57717E'));c.drawString(48,29,'SMRL | Zymo saved-run review | 2026-10-06');c.drawRightString(547,29,str(d.page))
SimpleDocTemplate(str(out),pagesize=(595,842),rightMargin=48,leftMargin=48,topMargin=42,bottomMargin=58,title='Zymo fecal experiment - TaxoViz and LRBinner',author='SMRL research workspace').build(story,onFirstPage=footer,onLaterPages=footer)
print(out);print('Context bytes:',cp.stat().st_size)