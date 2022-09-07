mkdir /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/TA_launch
cd /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/TA_launch

use UGER
use R-3.5

/xchip/beroukhimlab/Alex/PCAWG/Pan/Cluster/LaunchqsubUGER_20210111 -b /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/20220415_germsoma_sublist.sh --launch -t 30 -m 8
