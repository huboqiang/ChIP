git add .gitignore
git commit -m "Ignore file."

git add	__init__.py
git commit -m "Initiated file."

git add	bin/
git commit -m "Executable files for running."
git add	module00_StatInfo.py
git commit -m "File for sample information."

git add	module01_mapping_from_raw.py
git commit -m "Module1, generate scripts for mapping."

git add	module02_call_peaks.py
git commit -m "Module2, generate scripts for peak calling (mouse)."

git add	module02_call_peaks_human.py
git commit -m "Module2, generate scripts for peak calling (human)."

git add	module03_geneDensity.py
git commit -m "Module3, calculate read density on genes/XXbp Tiles."

git add	module04_ChIP_RNA_cor.py
git commit -m "Module4, RNA ChIP correlation on Promoter/Genebody."

git add	module_LoadSamp.py
git commit -m "File for sample information."

git add	module_Matrix.py
git commit -m "File for Matrix."

git add	module_PlotChIP_RNA.py
git commit -m "File for Plotting, 1."

git add	module_PlotRNA_RankedChip_gene.py
git commit -m "File for Plotting, 2."

git add	module_StatFunc.py
git commit -m "Submodule for statistics."
git add	module_refGene.py
git commit -m "Submodule for refGene.txt."
git add	module_running_jobs.py
git commit -m "Submodule for job running."









git add bin/PlotChip_IGV.py
git commit -m "plot chip-seq IGV like data."

git add module00_StatInfo.py 
git commit -m "StatInfo for samples"
git add module01_mapping_from_raw.py
git commit -m "modules for mapping pipeline."
git add module02_call_peaks.py
git commit -m "modules for calling peaks."
git add module03_geneDensity.py 
git commit -m "modules for gene density genome-wide or TSS region."

git rm module02_call_peaks_human.py 
git commit -m "remove module02_call_peaks_human.py previous version."

git rm module_LoadSamp.py 
git commit -m "remove module_LoadSamp.py       previous version."
git rm module_Matrix.py
git commit -m "remove module_Matrix.py         previous version."
git rm module_PlotChIP_RNA.py
git commit -m "remove module_PlotChIP_RNA.py   previous version."
git rm module_PlotRNA_RankedChip_gene.py
git commit -m "remove module_PlotRNA_RankedChip_gene.py      previous version."
git rm module_StatFunc.py
git commit -m "remove module_StatFunc.py        previous version."
git rm module_refGene.py
git commit -m "remove module_refGene.py         previous version."
git rm module_running_jobs.py
git commit -m "remove module_running_jobs.py    previous version."




git add	bin/div_bins/
git commit -m "calculate read signal density from BED file input."

git add	bin/find_ExonIntronIntergenic/
git commit -m "given a refGene file, define exon/intron/intergenic region for a genome."
git add	bin/multi-process.pl
git commit -m "perl scripts for running multiple shells in command line."
git add	bin/qsub-sge.pl
git commit -m "perl scripts for running multiple shells in SGE-qsub system."
git add	bin/s03_genePred2bed.py
git commit -m "convert gene-pred format into bed format."
git add	settings/
git commit -m "setting information for project-path and scripts for each step."
git add	utils/
git commit -m "modules for support."


git add run_chipseq.py
git commit -m "running the pipeline."