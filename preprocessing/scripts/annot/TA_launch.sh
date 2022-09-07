mkdir /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/annot/TA_launch
cd /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/annot/TA_launch

use UGER
use R-3.5

/xchip/beroukhimlab/Alex/PCAWG/Pan/Cluster/LaunchqsubUGER_20210111 -b /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/annot/20220419_annot_TA.sh --launch -t 30 -m 8
