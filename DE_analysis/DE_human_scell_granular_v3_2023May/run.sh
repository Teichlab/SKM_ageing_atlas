# working directory /nfs/team205/vk8/scripts/de_lmm_natsuhiko/DE_skm_scell_2021Jul
MEM=150000
if [ -z $1 ]
then
	bsub -q yesterday -e stderr -o stdout -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -n6 -G teichlab "sh run.sh farm"
	exit
fi
/software/team170/miniconda3/envs/seurat/bin/R --vanilla --quiet --args < run.R


#bsub -q long -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -Is $SHELL
