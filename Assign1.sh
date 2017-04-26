#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1


module load gcc openmpi
mpirun -np 1 ./assign1 1 100 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 4

module load gcc openmpi
mpirun -np 4 ./assign1 4 100 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -p node -n 25

module load gcc openmpi
mpirun -np 25 ./assign1 25 100 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1

module load gcc openmpi
mpirun -np 1 ./assign1 1 180 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 4

module load gcc openmpi
mpirun -np 4 ./assign1 4 180 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 9

module load gcc openmpi
mpirun -np 9 ./assign1 9 180 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -p node -n 36

module load gcc openmpi
mpirun -np 36 ./assign1 36 180 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1

module load gcc openmpi
mpirun -np 1 ./assign1 1 360 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 4

module load gcc openmpi
mpirun -np 4 ./assign1 4 360 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 9

module load gcc openmpi
mpirun -np 9 ./assign1 9 360 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -p node -n 36

module load gcc openmpi
mpirun -np 36 ./assign1 36 360 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1

module load gcc openmpi
mpirun -np 1 ./assign1 1 720 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 4

module load gcc openmpi
mpirun -np 4 ./assign1 4 720 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 9

module load gcc openmpi
mpirun -np 9 ./assign1 9 720 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -p node -n 36

module load gcc openmpi
mpirun -np 36 ./assign1 36 720 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1

module load gcc openmpi
mpirun -np 1 ./assign1 1 1440 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 4

module load gcc openmpi
mpirun -np 4 ./assign1 4 1440 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 9

module load gcc openmpi
mpirun -np 9 ./assign1 9 1440 

#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -p node -n 36

module load gcc openmpi
mpirun -np 36 ./assign1 36 1440 

