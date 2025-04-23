srun --partition=defq --nodelist=node006 --cpus-per-task=10 --mem=100G --exclusive --pty bash
module load singularity
#skipable mkdir -p $HOME/rstudio/tmp/var/run
singularity exec \
  --bind $HOME/rstudio/tmp:/tmp \
  --bind $HOME/rstudio/tmp/var:/var/lib/rstudio-server \
  --bind $HOME/rstudio/tmp/var/run:/var/run/rstudio-server \
  docker://bioconductor/bioconductor_docker:latest \
  rserver --auth-none=1 --auth-pam-helper-path=pam-helper --www-port=8787 --server-user=$USER


  #in_new_terminal
  ssh -N -L 8787:localhost:8787 -J mafechkar@ares.vumc.nl mafechkar@node006.cluster

  #open 
  http://localhost:8787
