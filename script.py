import os
from Bio import SeqIO
import pandas

rerun = False

if rerun:
    "python mOTUlizer/bin/mOTUlize.py -k /home/moritz/projects/0026_Chlorobi/data/checkms.tsv --fnas genomes.txt  --txt -c 20 --similarities chlorobis.tsv --MC 70 --Mc 5 --SC 50 --Sc 10  --prefix chlorobi_mOTU_ > chlorobi_mOTU.json"

with open("chlorobi_mOTU.json") as handle:
    motus = json.load(handle)

ref = []
newctg_id2oldctg_id = {}
for k in tqdm(os.listdir("data/AGs")):
    gen = list(SeqIO.parse("data/AGs/" + k + "/" + k + ".fna", "fasta"))
    for i,g in enumerate(gen):
        newctg_id2oldctg_id[g.id] = k + ":" + str(i).zfill(4)
        g.id = newctg_id2oldctg_id[g.id]
        g.description = ""
    ref += gen

SeqIO.write(ref, "data/reference.fna", "fasta")

if rerun:
    run("bowtie2-build --threads 20 data/reference.fna data/reference", shell = True)

read_fold = "/home/moritz/data/data_submit/reads/1Msubsample_mg_hc"

samples = {v.replace(".fastq.gz","").replace("1Msubsample_", "").split("_")[0] for v in os.listdir(read_fold)}

for s in tqdm(samples):
    mapper = """bowtie2 -x data/reference -1 {fold}/1Msubsample_{sample}_R1.fastq.gz -2 {fold}/1Msubsample_{sample}_R2.fastq.gz  -p20 -S data/maps/{sample}.sam  2> /dev/null
samtools view -u  data/maps/{sample}.sam | samtools sort -@ 20 > data/maps/{sample}.bam 2> /dev/null
samtools index data/maps/{sample}.bam  2> /dev/null
samtools flagstats data/maps/{sample}.bam > data/maps/{sample}.stats2> /dev/null
rm data/maps/{sample}.sam
""".format(sample = s, fold = read_fold)
    call(mapper, shell = True)

read_fold2 = "data/reads/"

samples2 = {v.replace(".fastq.gz","").replace("1Msubsample_", "").split("_")[0] for v in os.listdir(read_fold2) if v.startswith("1Msubsample_")}

for s in tqdm(samples2):
    mapper = """bowtie2 -x data/reference -1 {fold}/1Msubsample_{sample}_R1.fastq -2 {fold}/1Msubsample_{sample}_R2.fastq  -p20 -S data/maps/{sample}.sam  2> /dev/null
samtools view -u  data/maps/{sample}.sam | samtools sort -@ 20 > data/maps/{sample}.bam 2> /dev/null
samtools index data/maps/{sample}.bam 2> /dev/null
samtools flagstats data/maps/{sample}.bam > data/maps/{sample}.stats > /dev/null
rm data/maps/{sample}.sam
""".format(sample = s, fold = read_fold2)
    call(mapper, shell = True)

read_fold3 = "data/reads/metagenomes-wisc"
samples3 = {v.replace(".fastq","").replace("1Msubsample_", "").split("_")[0] for v in os.listdir(read_fold3) if v.startswith("1Msubsample_")}

for s in tqdm(samples3):
    mapper = """bowtie2 -x data/reference -1 {fold}/1Msubsample_{sample}.fastq -2 {fold}/1Msubsample_{sample}.fastq -p20 -S data/maps/{sample}.sam 2> /dev/null
samtools view -u  data/maps/{sample}.sam | samtools sort -@ 20 > data/maps/{sample}.bam 2> /dev/null
samtools index data/maps/{sample}.bam 2> /dev/null
samtools flagstats data/maps/{sample}.bam > data/maps/{sample}.stats 2> /dev/null
rm data/maps/{sample}.sam
""".format(sample = s, fold = read_fold3)
    call(mapper, shell = True)

"jgi_summarize_bam_contig_depths --percentIdentity 80 --noIntraDepthVariance  --referenceFasta data/reference.fna --outputDepth  data/cov_table.tsv  data/maps/*.bam"


covs = pandas.read_csv("data/cov_table.tsv", sep="\t", index_col=0)
del covs['totalAvgDepth']
lens = covs.contigLen.to_dict()
del covs['contigLen']
covs.columns = [c[:-4] for c in covs.columns]
tt = covs.transpose()
bins = {b.split(":")[0] for b in covs.index}

tt = { k : {kk : vv*lens[kk] for kk, vv in v.items()} for k, v in tqdm(tt.iterrows())}

sample2bin2reads = {s : {b : 0 for b in bins} for s in tt}
for k,v in tqdm(tt.items()):
    for kk, vv in v.items():
        sample2bin2reads[k][kk.split(":")[0]] += vv
bases = pandas.DataFrame.from_dict(sample2bin2reads)

libs_len = {'IHWH' : 218276636,
            'IHUI' : 201730774,
            'IHXT' : 223220788
            }

def nb_reads(s) :
    if s.startswith("SRR"):
        return 402000000
    if s.startswith("IH"):
        return libs_len[s]
    return 250000000

fraction = bases.apply(lambda x : pandas.core.series.Series([v/nb_reads(s) for v,s in zip(x, x.index)], index=x.index), axis=1)

motu2bin = {k : [vv['name'] for vv in v['genomes']] for k, v in motus.items()}
motu2bin['others'] = list(set(bins).difference(sum(list(motu2bin.values()),[])))
bin2motu = {vv : k for k,v in motu2bin.items() for vv in v }

sample2mOTU2fract = {k : {b : 0 for b in motu2bin} for k in fraction.columns}
for k, v in tqdm(fraction.iterrows()):
    for kk, vv in v.items():
        sample2mOTU2fract[kk][bin2motu.get(k, "other")] += vv

motu_wise = pandas.DataFrame.from_dict(sample2mOTU2fract)
motu_wise.to_csv("data/motu_table.csv")
out_groups = {k for k, v in motu2bin.items() if any([vv.startswith("GC") for vv in v])}
#cutoffs = motu_wise.loc[out_groups].apply(max).to_dict()
cutoff = max(fraction[[c for c in fraction.columns if c.startswith("SA2")]].sum())

motu_wise.apply(lambda x : pandas.core.series.Series([v if v > cutoffs[s] else 0.0 for v,s in zip(x, x.index)], index=x.index), axis=1)
motu_wise_cutoff = motu_wise.apply(lambda x : pandas.core.series.Series([v if v > cutoff else 0.0 for v,s in zip(x, x.index)], index=x.index), axis=1)
#motu_wise_cutoff = motu_wise.apply(lambda x : pandas.core.series.Series([v if v > cutoffs[s] else 0.0 for v,s in zip(x, x.index)], index=x.index), axis=1)
motu_wise.to_csv("data/motu_table.csv")
motu_wise_cutoff.to_csv("data/motu_table_cutoff.csv")

from Bio import Phylo
from ete3 import Tree
tree = list(Phylo.parse('data/gtdbtk_tree/gtdbtk_c__Chlorobia.tree', 'newick'))[0]

for c in tree.find_clades():
    if not c.is_terminal():
        c.name = None

Phylo.write(tree, "data/gtdbtk_tree/gtdbtk.tree", "newick")
ete_tree = Tree("data/gtdbtk_tree/gtdbtk.tree")

for l in ete_tree.iter_leaves():
    l.name = str(bin2motu.get(l.name, bin2motu.get(l.name[3:]))) +  " --- " + l.name
ete_tree.write(outfile = "data/gtdbtk_tree/gtdbtk_pretty.tree", "newick")

for f,v  in motu2bin.items():
    for vv in v:
        shutil.copy("data/faas/" + vv + ".faa", "data/mOTUs/" + f + "/" + vv + ".faa")

"python mOTUlizer/bin/mOTUpan.py --txt --faas <(find /home/moritz/projects/0026_Chlorobi/data/faas -name "*.faa") --genome2cog_only --checkm /home/moritz/projects/0026_Chlorobi/data/checkms.tsv -o all_chbi.json"

bin2oxid = {g : set() for g in bin2motu}
gene2genome = {v.id : f[:-4] for f in os.listdir("data/faas/") for v in SeqIO.parse("data/faas/" + f, "fasta")}


bin2oxid = {g : set() for g in bin2motu}
for f in os.listdir("data/gene-trees/"):
    gene = f[:-4]
    for v in SeqIO.parse("data/gene-trees/" + f, "fasta"):
        bin_id = gene2genome[v.id.replace("chlorobi-", "")]
        bin2oxid[bin_id].update({gene})

with open("data/genome2oxid.json", "w") as handle:
    json.dump({k : list(v) for k,v in bin2oxid.items()}, handle)

"for f in `ls /home/moritz/projects/0026_Chlorobi/data/mOTUs/`; do python mOTUlizer/bin/mOTUpan.py --txt --faas <(find /home/moritz/projects/0026_Chlorobi/data/mOTUs/$f/ -name "*.faa")  --checkm /home/moritz/projects/0026_Chlorobi/data/checkms.tsv -o /home/moritz/projects/0026_Chlorobi/data/dirtymotupan/$f.json --max_iter 1 --name $f --cog_file /home/moritz/projects/0026_Chlorobi/data/genome2oxid.json; done"

motu2oxid = dict()
for f in os.listdir("/home/moritz/projects/0026_Chlorobi/data/dirtymotupan/"):
    with open("/home/moritz/projects/0026_Chlorobi/data/dirtymotupan/" + f) as handle:
        dd = json.load(handle)
        motu2oxid[f[:-5]] = dd[f[:-5]]['likelies']

genes = pandas.DataFrame.from_dict(motu2oxid).fillna(-10).transpose()
genes['nb_genomes'] = [len(motu2bin[i]) for i in genes.index]
genes['genomes'] = [";".join(motu2bin[i]) for i in genes.index]
genes['mean_completeness'] = [ numpy.mean([checkms[vv]['Completeness'] for vv in  motu2bin[i]]) for i in genes.index]
genes['mean_contamination'] = [ numpy.mean([checkms[vv]['Contamination'] for vv in  motu2bin[i]]) for i in genes.index]

def get_genome_size(f):
    return sum([len(f.seq) for f in SeqIO.parse("data/fnas/" + f + ".fna", "fasta")])

genes['estimated_genome_size'] = [ numpy.mean([(100*get_genome_size(vv)/checkms[vv]['Completeness']) - (get_genome_size(vv)*checkms[vv]['Contamination']/100)  for vv in  motu2bin[i]]) for i in genes.index]
genes['mean_assembly_size'] = [ numpy.mean([get_genome_size(vv) for vv in  motu2bin[i]]) for i in genes.index]

bin_md = pandas.read_csv("/home/moritz/data/data_submit/metadata/master_table.csv", index_col=0)
bin_md.index = [b.replace("_megahit_metabat_", "_") for b in bin_md.index]
bin_md = bin_md.loc[[l for l in bin_md.index if l in bins]]
to_fill = {b : {} for b in bins.difference(bin_md.index)}
sub_bin_md = bin_md[['completeness','contamination','length']]
sub_bin_md['accession'] = None

bin2gca = {"Tsuji_L227_2013_bin56": "GCA_005843835.1",
           "Tsuji_L227_2013_bin55": "GCA_005843805.1",
           "Tsuji_L442_2011_2014_bin64" : "GCA_005843815.1",
           "Tsuji_L227_2013_bin22": "GCA_005843825.1",
           "Tsuji_L304_enr_S_6D_bin01": "GCA_005862285.1",
           "Tsuji_L227_enr_S_6D_bin01": "GCA_005862225.1",
           "TBepi_metabat_211": "IMG:2582580622",
           "TBhypo_metabat_111": "IMG:2582580651",
           "TBepi_metabat_2493": "IMG:2582580625",
           "TBhypo_metabat_03520": "IMG:2593339184"
           }

for bin, md in to_fill.items():
    seq = [s.seq for s in SeqIO.parse("data/fnas/" + bin + ".fna", "fasta")]
    md['nb_contigs'] = len(seq)
    md['completeness'] = checkms[bin]['Completeness']
    md['contamination'] = checkms[bin]['Contamination']
    md['length'] = sum([len(s) for s in seq])
    md['genome_accession'] = bin if bin.startswith("GC") else bin2gca.get(bin,"NA")
    md['GC'] = sum([s.count("G")+s.count("C") for s in seq])/md['length']
    seq = [s.seq for s in SeqIO.parse("data/faas/" + bin + ".faa", "fasta")]
    md['nb_proteins'] = len(seq)
    md['coding_density'] = 3*sum([len(s) for s in seq])/md['length']

to_fill.update(sub_bin_md.to_dict(orient = "index"))

for bin, md in to_fill.items():
    md['mOTU'] = bin2motu[bin]

sample_md = pandas.read_csv("/home/moritz/data/data_submit/metadata/samples_contextual_final.csv", index_col=0)
sample_md.index = [l.replace("_","-") for l in sample_md.index]

trouts = {"IHWH" : "IMG:3300020735",
"IHXT" : "IMG:3300020721",
"IHUI" : "IMG:3300020685"}

trout_date = {"IHWH" : "2007-07-25",
"IHXT" : "2008-07-22",
"IHUI" : "2009-08-11"}


kenora_md = { 'SRR8517163' : { 'date' : '2013-06-25', 'depth' : '6'},
              'SRR8517160' : { 'date' : '2014-07-08', 'depth' : '8'},
              'SRR8517158' : { 'date' : '2014-07-09', 'depth' : '13'},
              'SRR8517161' : { 'date' : '2014-07-08', 'depth' : '6'},
              'SRR8517159' : { 'date' : '2011-06-27', 'depth' : '16.5'},
              'SRR8517162' : { 'date' : '2013-06-25', 'depth' : '8'}
}

get_md = lambda s,md :  sample_md.loc[s if s[1] != '-' else s.replace("-",""), md]
all_samples = {s : {'dataset' : 'stratfreshdb', 'accession' : "SRA:" + get_md(s,"sample_accession"), 'date' : get_md(s,"collection date"), 'depth' : get_md(s, 'geographic location (depth)')} for s in samples if not s.startswith("SA") and not s.startswith("Mycocos")  and not s.startswith("SB-")}
all_samples.update( {s : {'dataset' : 'kenora', 'accession' : "SRA:" +  s, 'date' : kenora_md[s]['date'], 'depth' : kenora_md[s]['depth']} for s in samples2})
all_samples.update( {s : {'dataset' : 'troutbog', 'accession' : trouts[s], 'date' : trout_date[s], 'depth' : "NA"} for s in samples3})
 pandas.DataFrame.from_dict(all_samples, orient="index").to_csv("data/sample_data.csv", index_label="sample_ID")
print("Add tax and accession")


ag2acc = sag_dat['genome_accession'].to_dict()
